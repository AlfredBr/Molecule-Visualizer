/*
 * MolVis - Metal Compute Shaders
 *
 * GPU-accelerated molecular rendering using Metal compute shaders.
 * Renders ball-and-stick molecular models with Phong shading.
 *
 * Ported from CUDA implementation (cuda_renderer.cu)
 */

#include <metal_stdlib>
using namespace metal;

// ============== Constants ==============

// CPK colors for atoms (R, G, B)
constant float3 atomColors[18] = {
    float3(0.95f, 0.95f, 0.95f),  // H - white
    float3(0.2f,  0.2f,  0.2f),   // C - dark gray
    float3(0.2f,  0.3f,  0.9f),   // N - blue
    float3(0.9f,  0.2f,  0.2f),   // O - red
    float3(1.0f,  0.5f,  0.0f),   // P - orange
    float3(0.9f,  0.8f,  0.2f),   // S - yellow
    float3(0.2f,  0.9f,  0.2f),   // Cl - green
    float3(0.6f,  0.1f,  0.1f),   // Br - dark red
    float3(0.5f,  0.9f,  0.5f),   // F - light green
    float3(0.5f,  0.1f,  0.5f),   // I - purple
    float3(0.7f,  0.5f,  0.9f),   // Na - metallic purple
    float3(0.85f, 0.75f, 0.55f),  // Si - tan
    float3(1.0f,  0.65f, 0.65f),  // B - salmon
    float3(0.88f, 0.4f,  0.2f),   // Fe - orange/brown
    float3(0.85f, 0.55f, 0.2f),   // Cu - copper
    float3(0.75f, 0.75f, 0.8f),   // Al - silver
    float3(0.6f,  0.6f,  0.65f),  // Ti - gray
    float3(0.85f, 0.85f, 0.88f),  // Pt - white
};

// ============== Data Structures ==============

struct Atom {
    float x, y, z;
    int type;
    float radius;
};

struct Bond {
    int atom1, atom2;
    int order;
};

struct RenderParams {
    int width;
    int height;
    int numAtoms;
    int numBonds;
    float rotX;
    float rotY;
    float zoom;
    float offsetX;
    float offsetY;
};

struct TextParams {
    int textLen;
    int startX;
    int startY;
    int scale;
    int width;
    int height;
};

// ============== 3D Math Helpers ==============

inline float3 normalize3(float3 v) {
    float len = length(v);
    if (len > 0.0001f) {
        return v / len;
    }
    return float3(0, 1, 0);
}

inline float dot3(float3 a, float3 b) {
    return dot(a, b);
}

inline float3 reflect3(float3 v, float3 n) {
    float d = 2.0f * dot(v, n);
    return v - d * n;
}

// Ray-sphere intersection
inline float intersectSphere(float3 ro, float3 rd, float3 center, float radius) {
    float3 oc = ro - center;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - radius * radius;
    float h = b * b - c;
    if (h < 0.0f) return -1.0f;
    return -b - sqrt(h);
}

// Ray-cylinder intersection
inline float intersectCylinder(float3 ro, float3 rd, float3 pa, float3 pb, float radius, thread float3* outNormal) {
    float3 ba = pb - pa;
    float3 oc = ro - pa;

    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baoc = dot(ba, oc);

    float k2 = baba - bard * bard;
    float k1 = baba * dot(oc, rd) - baoc * bard;
    float k0 = baba * dot(oc, oc) - baoc * baoc - radius * radius * baba;

    if (abs(k2) < 0.0001f) return -1.0f;

    float h = k1 * k1 - k2 * k0;
    if (h < 0.0f) return -1.0f;

    h = sqrt(h);
    float t = (-k1 - h) / k2;

    float y = baoc + t * bard;
    if (y > 0.0f && y < baba && t > 0.0f) {
        float3 hitPoint = ro + t * rd;
        float3 toHit = hitPoint - pa;
        float proj = dot(toHit, ba) / baba;
        float3 onAxis = pa + proj * ba;
        *outNormal = normalize(hitPoint - onAxis);
        return t;
    }

    return -1.0f;
}

// ============== Clear Kernel ==============

kernel void clearKernel(
    texture2d<half, access::write> output [[texture(0)]],
    uint2 gid [[thread_position_in_grid]])
{
    uint width = output.get_width();
    uint height = output.get_height();

    if (gid.x >= width || gid.y >= height) return;

    // Gradient background (dark blue-gray)
    float gy = float(gid.y) / float(height);
    half4 color = half4(
        half((20.0f + gy * 30.0f) / 255.0f),
        half((25.0f + gy * 35.0f) / 255.0f),
        half((35.0f + gy * 40.0f) / 255.0f),
        1.0h
    );
    output.write(color, gid);
}

// ============== Bitmap Font for Text Rendering ==============

// 6x8 bitmap font (ASCII 32-127) - stored as constant arrays
constant uchar fontData[96][8] = {
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // 32 space
    {0x04,0x04,0x04,0x04,0x04,0x00,0x04,0x00}, // 33 !
    {0x0A,0x0A,0x0A,0x00,0x00,0x00,0x00,0x00}, // 34 "
    {0x0A,0x0A,0x1F,0x0A,0x1F,0x0A,0x0A,0x00}, // 35 #
    {0x04,0x0F,0x14,0x0E,0x05,0x1E,0x04,0x00}, // 36 $
    {0x18,0x19,0x02,0x04,0x08,0x13,0x03,0x00}, // 37 %
    {0x08,0x14,0x14,0x08,0x15,0x12,0x0D,0x00}, // 38 &
    {0x04,0x04,0x08,0x00,0x00,0x00,0x00,0x00}, // 39 '
    {0x02,0x04,0x08,0x08,0x08,0x04,0x02,0x00}, // 40 (
    {0x08,0x04,0x02,0x02,0x02,0x04,0x08,0x00}, // 41 )
    {0x00,0x04,0x15,0x0E,0x15,0x04,0x00,0x00}, // 42 *
    {0x00,0x04,0x04,0x1F,0x04,0x04,0x00,0x00}, // 43 +
    {0x00,0x00,0x00,0x00,0x00,0x04,0x04,0x08}, // 44 ,
    {0x00,0x00,0x00,0x1F,0x00,0x00,0x00,0x00}, // 45 -
    {0x00,0x00,0x00,0x00,0x00,0x00,0x04,0x00}, // 46 .
    {0x01,0x02,0x02,0x04,0x08,0x08,0x10,0x00}, // 47 /
    {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E,0x00}, // 48 0
    {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E,0x00}, // 49 1
    {0x0E,0x11,0x01,0x02,0x04,0x08,0x1F,0x00}, // 50 2
    {0x0E,0x11,0x01,0x06,0x01,0x11,0x0E,0x00}, // 51 3
    {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02,0x00}, // 52 4
    {0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E,0x00}, // 53 5
    {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E,0x00}, // 54 6
    {0x1F,0x01,0x02,0x04,0x08,0x08,0x08,0x00}, // 55 7
    {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E,0x00}, // 56 8
    {0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C,0x00}, // 57 9
    {0x00,0x00,0x04,0x00,0x00,0x04,0x00,0x00}, // 58 :
    {0x00,0x00,0x04,0x00,0x00,0x04,0x04,0x08}, // 59 ;
    {0x02,0x04,0x08,0x10,0x08,0x04,0x02,0x00}, // 60 <
    {0x00,0x00,0x1F,0x00,0x1F,0x00,0x00,0x00}, // 61 =
    {0x08,0x04,0x02,0x01,0x02,0x04,0x08,0x00}, // 62 >
    {0x0E,0x11,0x01,0x02,0x04,0x00,0x04,0x00}, // 63 ?
    {0x0E,0x11,0x17,0x15,0x17,0x10,0x0E,0x00}, // 64 @
    {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11,0x00}, // 65 A
    {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E,0x00}, // 66 B
    {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E,0x00}, // 67 C
    {0x1E,0x11,0x11,0x11,0x11,0x11,0x1E,0x00}, // 68 D
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F,0x00}, // 69 E
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10,0x00}, // 70 F
    {0x0E,0x11,0x10,0x17,0x11,0x11,0x0F,0x00}, // 71 G
    {0x11,0x11,0x11,0x1F,0x11,0x11,0x11,0x00}, // 72 H
    {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E,0x00}, // 73 I
    {0x07,0x02,0x02,0x02,0x02,0x12,0x0C,0x00}, // 74 J
    {0x11,0x12,0x14,0x18,0x14,0x12,0x11,0x00}, // 75 K
    {0x10,0x10,0x10,0x10,0x10,0x10,0x1F,0x00}, // 76 L
    {0x11,0x1B,0x15,0x15,0x11,0x11,0x11,0x00}, // 77 M
    {0x11,0x19,0x15,0x13,0x11,0x11,0x11,0x00}, // 78 N
    {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E,0x00}, // 79 O
    {0x1E,0x11,0x11,0x1E,0x10,0x10,0x10,0x00}, // 80 P
    {0x0E,0x11,0x11,0x11,0x15,0x12,0x0D,0x00}, // 81 Q
    {0x1E,0x11,0x11,0x1E,0x14,0x12,0x11,0x00}, // 82 R
    {0x0E,0x11,0x10,0x0E,0x01,0x11,0x0E,0x00}, // 83 S
    {0x1F,0x04,0x04,0x04,0x04,0x04,0x04,0x00}, // 84 T
    {0x11,0x11,0x11,0x11,0x11,0x11,0x0E,0x00}, // 85 U
    {0x11,0x11,0x11,0x11,0x11,0x0A,0x04,0x00}, // 86 V
    {0x11,0x11,0x11,0x15,0x15,0x1B,0x11,0x00}, // 87 W
    {0x11,0x11,0x0A,0x04,0x0A,0x11,0x11,0x00}, // 88 X
    {0x11,0x11,0x0A,0x04,0x04,0x04,0x04,0x00}, // 89 Y
    {0x1F,0x01,0x02,0x04,0x08,0x10,0x1F,0x00}, // 90 Z
    {0x0E,0x08,0x08,0x08,0x08,0x08,0x0E,0x00}, // 91 [
    {0x10,0x08,0x08,0x04,0x02,0x02,0x01,0x00}, // 92 backslash
    {0x0E,0x02,0x02,0x02,0x02,0x02,0x0E,0x00}, // 93 ]
    {0x04,0x0A,0x11,0x00,0x00,0x00,0x00,0x00}, // 94 ^
    {0x00,0x00,0x00,0x00,0x00,0x00,0x1F,0x00}, // 95 _
    {0x08,0x04,0x02,0x00,0x00,0x00,0x00,0x00}, // 96 `
    {0x00,0x00,0x0E,0x01,0x0F,0x11,0x0F,0x00}, // 97 a
    {0x10,0x10,0x1E,0x11,0x11,0x11,0x1E,0x00}, // 98 b
    {0x00,0x00,0x0E,0x11,0x10,0x11,0x0E,0x00}, // 99 c
    {0x01,0x01,0x0F,0x11,0x11,0x11,0x0F,0x00}, // 100 d
    {0x00,0x00,0x0E,0x11,0x1F,0x10,0x0E,0x00}, // 101 e
    {0x02,0x05,0x04,0x0E,0x04,0x04,0x04,0x00}, // 102 f
    {0x00,0x00,0x0F,0x11,0x11,0x0F,0x01,0x0E}, // 103 g
    {0x10,0x10,0x16,0x19,0x11,0x11,0x11,0x00}, // 104 h
    {0x04,0x00,0x0C,0x04,0x04,0x04,0x0E,0x00}, // 105 i
    {0x02,0x00,0x06,0x02,0x02,0x02,0x12,0x0C}, // 106 j
    {0x10,0x10,0x12,0x14,0x18,0x14,0x12,0x00}, // 107 k
    {0x0C,0x04,0x04,0x04,0x04,0x04,0x0E,0x00}, // 108 l
    {0x00,0x00,0x1A,0x15,0x15,0x11,0x11,0x00}, // 109 m
    {0x00,0x00,0x16,0x19,0x11,0x11,0x11,0x00}, // 110 n
    {0x00,0x00,0x0E,0x11,0x11,0x11,0x0E,0x00}, // 111 o
    {0x00,0x00,0x1E,0x11,0x11,0x1E,0x10,0x10}, // 112 p
    {0x00,0x00,0x0F,0x11,0x11,0x0F,0x01,0x01}, // 113 q
    {0x00,0x00,0x16,0x19,0x10,0x10,0x10,0x00}, // 114 r
    {0x00,0x00,0x0E,0x10,0x0E,0x01,0x1E,0x00}, // 115 s
    {0x04,0x04,0x0E,0x04,0x04,0x05,0x02,0x00}, // 116 t
    {0x00,0x00,0x11,0x11,0x11,0x13,0x0D,0x00}, // 117 u
    {0x00,0x00,0x11,0x11,0x11,0x0A,0x04,0x00}, // 118 v
    {0x00,0x00,0x11,0x11,0x15,0x15,0x0A,0x00}, // 119 w
    {0x00,0x00,0x11,0x0A,0x04,0x0A,0x11,0x00}, // 120 x
    {0x00,0x00,0x11,0x11,0x11,0x0F,0x01,0x0E}, // 121 y
    {0x00,0x00,0x1F,0x02,0x04,0x08,0x1F,0x00}, // 122 z
    {0x02,0x04,0x04,0x08,0x04,0x04,0x02,0x00}, // 123 {
    {0x04,0x04,0x04,0x04,0x04,0x04,0x04,0x00}, // 124 |
    {0x08,0x04,0x04,0x02,0x04,0x04,0x08,0x00}, // 125 }
    {0x00,0x00,0x08,0x15,0x02,0x00,0x00,0x00}, // 126 ~
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // 127
};

inline uchar getPixel(int c, int px, int py) {
    if (c < 32 || c > 127) return 0;
    int idx = c - 32;
    if (idx < 0 || idx >= 96) return 0;
    if (py < 0 || py >= 8) return 0;
    if (px < 0 || px >= 6) return 0;

    uchar row = fontData[idx][py];
    return (row >> (4 - px)) & 1;
}

// Helper functions for subscript detection (chemical formulas)
inline bool isElementChar(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

inline bool isDigit(char c) {
    return c >= '0' && c <= '9';
}

inline bool shouldSubscript(device const char* text, int textLen, int charIndex) {
    if (charIndex <= 0 || charIndex >= textLen) return false;
    char current = text[charIndex];
    char prev = text[charIndex - 1];

    if (isDigit(current) && isElementChar(prev)) return true;

    if (isDigit(current) && isDigit(prev) && charIndex >= 2) {
        int i = charIndex - 1;
        while (i >= 0 && isDigit(text[i])) i--;
        if (i >= 0 && isElementChar(text[i])) return true;
    }

    return false;
}

// ============== Text Rendering Kernel ==============

kernel void renderTextKernel(
    texture2d<half, access::read_write> output [[texture(0)]],
    device const char* text [[buffer(0)]],
    constant TextParams& params [[buffer(1)]],
    uint2 gid [[thread_position_in_grid]])
{
    int x = gid.x;
    int y = gid.y;

    if (x >= params.width || y >= params.height) return;

    int charHeight = 8 * params.scale;
    int subScale = (params.scale * 2) / 3;
    if (subScale < 1) subScale = 1;
    int subCharWidth = 6 * subScale;
    int subCharHeight = 8 * subScale;
    int subOffsetY = charHeight - subCharHeight;

    int curX = params.startX;
    int hitChar = -1;
    int charStartX = 0;
    bool hitIsSubscript = false;

    for (int i = 0; i < params.textLen; i++) {
        bool isSub = shouldSubscript(text, params.textLen, i);
        int cw = isSub ? subCharWidth : (6 * params.scale);

        if (x >= curX && x < curX + cw) {
            hitChar = i;
            charStartX = curX;
            hitIsSubscript = isSub;
            break;
        }
        curX += cw;
    }

    if (hitChar < 0) return;

    int charY, charH, charScale;
    if (hitIsSubscript) {
        charY = params.startY + subOffsetY;
        charH = subCharHeight;
        charScale = subScale;
    } else {
        charY = params.startY;
        charH = charHeight;
        charScale = params.scale;
    }

    if (y < charY || y >= charY + charH) return;

    int relX = x - charStartX;
    int relY = y - charY;
    int pixelX = relX / charScale;
    int pixelY = relY / charScale;

    uchar c = text[hitChar];
    if (getPixel(c, pixelX, pixelY)) {
        // White text
        output.write(half4(1.0h, 1.0h, 1.0h, 1.0h), gid);
    }
}

// ============== Molecule Rendering Kernel ==============

kernel void renderMoleculeKernel(
    texture2d<half, access::read_write> output [[texture(0)]],
    device const Atom* atoms [[buffer(0)]],
    device const Bond* bonds [[buffer(1)]],
    constant RenderParams& params [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]])
{
    int x = gid.x;
    int y = gid.y;

    if (x >= params.width || y >= params.height) return;

    // Camera setup
    float aspectRatio = float(params.width) / float(params.height);
    float fovScale = tan(0.5f * 0.8f);

    float ndcX = (2.0f * float(x) / float(params.width) - 1.0f) * aspectRatio * fovScale - params.offsetX;
    float ndcY = (1.0f - 2.0f * float(y) / float(params.height)) * fovScale - params.offsetY;

    // Ray direction
    float3 rd = normalize(float3(ndcX, ndcY, 1.0f));

    // Rotate ray
    float cosY = cos(params.rotY), sinY = sin(params.rotY);
    float cosX = cos(params.rotX), sinX = sin(params.rotX);

    float rx = rd.x * cosY + rd.z * sinY;
    float rz = -rd.x * sinY + rd.z * cosY;
    rd.x = rx; rd.z = rz;

    float ry = rd.y * cosX - rd.z * sinX;
    rz = rd.y * sinX + rd.z * cosX;
    rd.y = ry; rd.z = rz;

    // Ray origin
    float3 ro = float3(0, 0, -params.zoom);

    // Rotate camera
    float cox = ro.x * cosY + ro.z * sinY;
    float coz = -ro.x * sinY + ro.z * cosY;
    ro.x = cox; ro.z = coz;

    float coy = ro.y * cosX - ro.z * sinX;
    coz = ro.y * sinX + ro.z * cosX;
    ro.y = coy; ro.z = coz;

    // Lighting
    float3 lightDir = normalize(float3(0.5f, 0.8f, -0.3f));
    float3 lightDir2 = normalize(float3(-0.3f, 0.2f, 0.5f));

    float minT = 1e10f;
    float3 hitNormal = float3(0, 0, 0);
    float3 hitColor = float3(0, 0, 0);
    int hitType = 0;  // 0=none, 1=atom, 2=bond

    // Test atoms (spheres)
    for (int i = 0; i < params.numAtoms; i++) {
        float3 center = float3(atoms[i].x, atoms[i].y, atoms[i].z);
        float radius = atoms[i].radius;

        float t = intersectSphere(ro, rd, center, radius);
        if (t > 0.0f && t < minT) {
            minT = t;
            float3 hitPos = ro + t * rd;
            hitNormal = normalize(hitPos - center);
            int atomType = atoms[i].type;
            if (atomType >= 0 && atomType < 18) {
                hitColor = atomColors[atomType];
            } else {
                hitColor = float3(0.5f, 0.5f, 0.5f);
            }
            hitType = 1;
        }
    }

    // Test bonds (cylinders)
    float bondRadius = 0.08f;
    for (int i = 0; i < params.numBonds; i++) {
        int a1 = bonds[i].atom1;
        int a2 = bonds[i].atom2;

        if (a1 < 0 || a1 >= params.numAtoms || a2 < 0 || a2 >= params.numAtoms) continue;

        float3 p1 = float3(atoms[a1].x, atoms[a1].y, atoms[a1].z);
        float3 p2 = float3(atoms[a2].x, atoms[a2].y, atoms[a2].z);

        int order = bonds[i].order;
        float3 bondDir = normalize(p2 - p1);

        // Find perpendicular vector
        float3 perp;
        if (abs(bondDir.x) < 0.9f) {
            perp = normalize(float3(0, -bondDir.z, bondDir.y));
        } else {
            perp = normalize(float3(-bondDir.z, 0, bondDir.x));
        }

        for (int b = 0; b < order; b++) {
            float3 offset = float3(0, 0, 0);
            if (order == 2) {
                float off = (b == 0) ? -0.08f : 0.08f;
                offset = perp * off;
            } else if (order == 3) {
                float off = float(b - 1) * 0.1f;
                offset = perp * off;
            }

            float3 bp1 = p1 + offset;
            float3 bp2 = p2 + offset;

            float3 cylNormal;
            float t = intersectCylinder(ro, rd, bp1, bp2, bondRadius, &cylNormal);
            if (t > 0.0f && t < minT) {
                minT = t;
                hitNormal = cylNormal;
                hitColor = float3(0.6f, 0.6f, 0.6f);
                hitType = 2;
            }
        }
    }

    if (hitType > 0) {
        // Phong shading
        float diffuse1 = max(0.0f, dot(hitNormal, lightDir));
        float diffuse2 = max(0.0f, dot(hitNormal, lightDir2)) * 0.3f;

        float3 viewDir = normalize(-rd);
        float3 reflectDir = reflect(-lightDir, hitNormal);
        float spec = pow(max(0.0f, dot(viewDir, reflectDir)), 40.0f);

        float ambient = 0.15f;
        float3 finalColor;

        if (hitType == 1) {
            // Atom: metallic look with environment reflection
            float3 envReflect = reflect(rd, hitNormal);
            float envBrightness = 0.5f + 0.5f * envReflect.y;

            finalColor = hitColor * (ambient + diffuse1 * 0.7f + diffuse2);
            finalColor += float3(0.8f, 0.85f, 0.9f) * spec * 0.4f;
            finalColor += hitColor * envBrightness * 0.15f;

            // Rim lighting
            float rim = 1.0f - max(0.0f, dot(viewDir, hitNormal));
            rim = pow(rim, 3.0f);
            finalColor += float3(0.3f, 0.35f, 0.4f) * rim * 0.3f;
        } else {
            // Bond: matte gray
            finalColor = hitColor * (ambient + diffuse1 * 0.8f + diffuse2);
            finalColor += float3(0.5f) * spec * 0.2f;
        }

        // Clamp and output
        finalColor = clamp(finalColor, 0.0f, 1.0f);
        output.write(half4(half3(finalColor), 1.0h), gid);
    }
}
