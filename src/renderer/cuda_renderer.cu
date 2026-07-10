/*
 * MolVis - CUDA Renderer Implementation
 *
 * GPU-accelerated molecular rendering with DX11-CUDA interop.
 * Renders ball-and-stick molecular models with Phong shading.
 */

#include "cuda_renderer.h"
#include <stdio.h>
#include <string.h>

#define PI 3.14159265359f

// ============== CUDA Device Constants ==============

// CPK colors for atoms (R, G, B)
__device__ __constant__ float3 atomColors[ATOM_TYPE_COUNT] = {
    {0.95f, 0.95f, 0.95f},  // H - white
    {0.2f,  0.2f,  0.2f},   // C - dark gray
    {0.2f,  0.3f,  0.9f},   // N - blue
    {0.9f,  0.2f,  0.2f},   // O - red
    {1.0f,  0.5f,  0.0f},   // P - orange
    {0.9f,  0.8f,  0.2f},   // S - yellow
    {0.2f,  0.9f,  0.2f},   // Cl - green
    {0.6f,  0.1f,  0.1f},   // Br - dark red
    {0.5f,  0.9f,  0.5f},   // F - light green
    {0.5f,  0.1f,  0.5f},   // I - purple
    {0.7f,  0.5f,  0.9f},   // Na - metallic purple
    {0.85f, 0.75f, 0.55f},  // Si - tan
    {1.0f,  0.65f, 0.65f},  // B - salmon
    {0.88f, 0.4f,  0.2f},   // Fe - orange/brown
    {0.85f, 0.55f, 0.2f},   // Cu - copper
    {0.75f, 0.75f, 0.8f},   // Al - silver
    {0.6f,  0.6f,  0.65f},  // Ti - gray
    {0.85f, 0.85f, 0.88f},  // Pt - white
    {0.51f, 0.51f, 0.56f},  // Re - metallic gray
    {0.26f, 0.62f, 0.69f},  // Xe - cyan (noble gas)
    {0.33f, 0.71f, 0.71f},  // Mo - teal
    {0.13f, 0.58f, 0.84f},  // W - blue
    {1.00f, 0.63f, 0.00f},  // Se - orange
    {0.40f, 0.56f, 0.56f},  // Ge - gray-green
    {0.76f, 0.56f, 0.56f},  // Ga - muted rose
};

// ============== 3D Math Helpers ==============

__device__ float3 normalize3(float3 v) {
    float len = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len > 0.0001f) {
        return make_float3(v.x / len, v.y / len, v.z / len);
    }
    return make_float3(0, 1, 0);
}

__device__ float dot3(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ float3 reflect3(float3 v, float3 n) {
    float d = 2.0f * dot3(v, n);
    return make_float3(v.x - d * n.x, v.y - d * n.y, v.z - d * n.z);
}

// Ray-sphere intersection
__device__ float intersectSphere(float3 ro, float3 rd, float3 center, float radius) {
    float3 oc = make_float3(ro.x - center.x, ro.y - center.y, ro.z - center.z);
    float b = dot3(oc, rd);
    float c = dot3(oc, oc) - radius * radius;
    float h = b * b - c;
    if (h < 0.0f) return -1.0f;
    return -b - sqrtf(h);
}

// Ray-cylinder intersection
__device__ float intersectCylinder(float3 ro, float3 rd, float3 pa, float3 pb, float radius, float3* outNormal) {
    float3 ba = make_float3(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
    float3 oc = make_float3(ro.x - pa.x, ro.y - pa.y, ro.z - pa.z);

    float baba = dot3(ba, ba);
    float bard = dot3(ba, rd);
    float baoc = dot3(ba, oc);

    float k2 = baba - bard * bard;
    float k1 = baba * dot3(oc, rd) - baoc * bard;
    float k0 = baba * dot3(oc, oc) - baoc * baoc - radius * radius * baba;

    if (fabsf(k2) < 0.0001f) return -1.0f;

    float h = k1 * k1 - k2 * k0;
    if (h < 0.0f) return -1.0f;

    h = sqrtf(h);
    float t = (-k1 - h) / k2;

    float y = baoc + t * bard;
    if (y > 0.0f && y < baba && t > 0.0f) {
        float3 hitPoint = make_float3(ro.x + t * rd.x, ro.y + t * rd.y, ro.z + t * rd.z);
        float3 toHit = make_float3(hitPoint.x - pa.x, hitPoint.y - pa.y, hitPoint.z - pa.z);
        float proj = dot3(toHit, ba) / baba;
        float3 onAxis = make_float3(pa.x + proj * ba.x, pa.y + proj * ba.y, pa.z + proj * ba.z);
        *outNormal = normalize3(make_float3(hitPoint.x - onAxis.x, hitPoint.y - onAxis.y, hitPoint.z - onAxis.z));
        return t;
    }

    return -1.0f;
}

// ============== CUDA Kernels ==============

__global__ void clearKernel(unsigned char* pixels, int width, int height) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) return;

    int idx = (y * width + x) * 4;

    // Gradient background (dark blue-gray)
    float gy = (float)y / height;
    pixels[idx + 0] = (unsigned char)(20 + gy * 30);   // R
    pixels[idx + 1] = (unsigned char)(25 + gy * 35);   // G
    pixels[idx + 2] = (unsigned char)(35 + gy * 40);   // B
    pixels[idx + 3] = 255;                              // A
}

// ============== Bitmap Font for Text Rendering ==============

// 6x8 bitmap font (ASCII 32-127)
__device__ unsigned char getPixel(int c, int px, int py) {
    if (c < 32 || c > 127) return 0;

    static const unsigned char font[][8] = {
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

    int idx = c - 32;
    if (idx < 0 || idx >= 96) return 0;
    if (py < 0 || py >= 8) return 0;
    if (px < 0 || px >= 6) return 0;

    unsigned char row = font[idx][py];
    return (row >> (4 - px)) & 1;
}

// Helper functions for subscript detection (chemical formulas)
__device__ bool isElementChar(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

__device__ bool isDigit(char c) {
    return c >= '0' && c <= '9';
}

__device__ bool shouldSubscript(const char* text, int textLen, int charIndex) {
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

// Text rendering kernel with subscript support for chemical formulas
__global__ void renderTextKernel(
    unsigned char* pixels, int width, int height,
    const char* text, int textLen, int startX, int startY, int scale)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) return;

    int charHeight = 8 * scale;
    int subScale = (scale * 2) / 3;
    if (subScale < 1) subScale = 1;
    int subCharWidth = 6 * subScale;
    int subCharHeight = 8 * subScale;
    int subOffsetY = charHeight - subCharHeight;

    int curX = startX;
    int hitChar = -1;
    int charStartX = 0;
    bool hitIsSubscript = false;

    for (int i = 0; i < textLen; i++) {
        bool isSub = shouldSubscript(text, textLen, i);
        int cw = isSub ? subCharWidth : (6 * scale);

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
        charY = startY + subOffsetY;
        charH = subCharHeight;
        charScale = subScale;
    } else {
        charY = startY;
        charH = charHeight;
        charScale = scale;
    }

    if (y < charY || y >= charY + charH) return;

    int relX = x - charStartX;
    int relY = y - charY;
    int pixelX = relX / charScale;
    int pixelY = relY / charScale;

    unsigned char c = text[hitChar];
    if (getPixel(c, pixelX, pixelY)) {
        int idx = (y * width + x) * 4;
        // White text
        pixels[idx + 0] = 255;  // R
        pixels[idx + 1] = 255;  // G
        pixels[idx + 2] = 255;  // B
        pixels[idx + 3] = 255;  // A
    }
}

__global__ void renderMoleculeKernel(
    unsigned char* pixels, int width, int height,
    Atom* atoms, int numAtoms,
    Bond* bonds, int numBonds,
    float rotX, float rotY, float zoom,
    float offsetX, float offsetY)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) return;

    // Camera setup
    float aspectRatio = (float)width / height;
    float fovScale = tanf(0.5f * 0.8f);

    float ndcX = (2.0f * x / width - 1.0f) * aspectRatio * fovScale - offsetX;
    float ndcY = (1.0f - 2.0f * y / height) * fovScale - offsetY;

    // Ray direction
    float3 rd = normalize3(make_float3(ndcX, ndcY, 1.0f));

    // Rotate ray
    float cosY = cosf(rotY), sinY = sinf(rotY);
    float cosX = cosf(rotX), sinX = sinf(rotX);

    float rx = rd.x * cosY + rd.z * sinY;
    float rz = -rd.x * sinY + rd.z * cosY;
    rd.x = rx; rd.z = rz;

    float ry = rd.y * cosX - rd.z * sinX;
    rz = rd.y * sinX + rd.z * cosX;
    rd.y = ry; rd.z = rz;

    // Ray origin
    float3 ro = make_float3(0, 0, -zoom);

    // Rotate camera
    float cox = ro.x * cosY + ro.z * sinY;
    float coz = -ro.x * sinY + ro.z * cosY;
    ro.x = cox; ro.z = coz;

    float coy = ro.y * cosX - ro.z * sinX;
    coz = ro.y * sinX + ro.z * cosX;
    ro.y = coy; ro.z = coz;

    // Lighting
    float3 lightDir = normalize3(make_float3(0.5f, 0.8f, -0.3f));
    float3 lightDir2 = normalize3(make_float3(-0.3f, 0.2f, 0.5f));

    float minT = 1e10f;
    float3 hitNormal;
    float3 hitColor;
    int hitType = 0;  // 0=none, 1=atom, 2=bond

    // Test atoms (spheres)
    for (int i = 0; i < numAtoms; i++) {
        float3 center = make_float3(atoms[i].x, atoms[i].y, atoms[i].z);
        float radius = atoms[i].radius;

        float t = intersectSphere(ro, rd, center, radius);
        if (t > 0.0f && t < minT) {
            minT = t;
            float3 hitPos = make_float3(ro.x + t * rd.x, ro.y + t * rd.y, ro.z + t * rd.z);
            hitNormal = normalize3(make_float3(hitPos.x - center.x, hitPos.y - center.y, hitPos.z - center.z));
            hitColor = atomColors[atoms[i].type];
            hitType = 1;
        }
    }

    // Test bonds (cylinders)
    float bondRadius = 0.08f;
    for (int i = 0; i < numBonds; i++) {
        int a1 = bonds[i].atom1;
        int a2 = bonds[i].atom2;

        float3 p1 = make_float3(atoms[a1].x, atoms[a1].y, atoms[a1].z);
        float3 p2 = make_float3(atoms[a2].x, atoms[a2].y, atoms[a2].z);

        int order = bonds[i].order;
        float3 bondDir = normalize3(make_float3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z));

        // Find perpendicular vector
        float3 perp;
        if (fabsf(bondDir.x) < 0.9f) {
            perp = normalize3(make_float3(0, -bondDir.z, bondDir.y));
        } else {
            perp = normalize3(make_float3(-bondDir.z, 0, bondDir.x));
        }

        for (int b = 0; b < order; b++) {
            float3 offset = make_float3(0, 0, 0);
            if (order == 2) {
                float off = (b == 0) ? -0.08f : 0.08f;
                offset = make_float3(perp.x * off, perp.y * off, perp.z * off);
            } else if (order == 3) {
                float off = (b - 1) * 0.1f;
                offset = make_float3(perp.x * off, perp.y * off, perp.z * off);
            }

            float3 bp1 = make_float3(p1.x + offset.x, p1.y + offset.y, p1.z + offset.z);
            float3 bp2 = make_float3(p2.x + offset.x, p2.y + offset.y, p2.z + offset.z);

            float3 cylNormal;
            float t = intersectCylinder(ro, rd, bp1, bp2, bondRadius, &cylNormal);
            if (t > 0.0f && t < minT) {
                minT = t;
                hitNormal = cylNormal;
                hitColor = make_float3(0.6f, 0.6f, 0.6f);
                hitType = 2;
            }
        }
    }

    if (hitType > 0) {
        // Phong shading
        float diffuse1 = fmaxf(0.0f, dot3(hitNormal, lightDir));
        float diffuse2 = fmaxf(0.0f, dot3(hitNormal, lightDir2)) * 0.3f;

        float3 viewDir = normalize3(make_float3(-rd.x, -rd.y, -rd.z));
        float3 reflectDir = reflect3(make_float3(-lightDir.x, -lightDir.y, -lightDir.z), hitNormal);
        float spec = powf(fmaxf(0.0f, dot3(viewDir, reflectDir)), 40.0f);

        float ambient = 0.15f;
        float3 finalColor;

        if (hitType == 1) {
            finalColor.x = hitColor.x * (ambient + diffuse1 * 0.7f + diffuse2) + spec * 0.6f;
            finalColor.y = hitColor.y * (ambient + diffuse1 * 0.7f + diffuse2) + spec * 0.6f;
            finalColor.z = hitColor.z * (ambient + diffuse1 * 0.7f + diffuse2) + spec * 0.6f;
        } else {
            finalColor.x = hitColor.x * (ambient + diffuse1 * 0.5f + diffuse2) + spec * 0.4f;
            finalColor.y = hitColor.y * (ambient + diffuse1 * 0.5f + diffuse2) + spec * 0.4f;
            finalColor.z = hitColor.z * (ambient + diffuse1 * 0.5f + diffuse2) + spec * 0.4f;
        }

        // Fresnel rim lighting
        float fresnel = powf(1.0f - fmaxf(0.0f, dot3(viewDir, hitNormal)), 3.0f);
        finalColor.x += fresnel * 0.15f;
        finalColor.y += fresnel * 0.15f;
        finalColor.z += fresnel * 0.2f;

        int idx = (y * width + x) * 4;
        pixels[idx + 0] = (unsigned char)fminf(255.0f, finalColor.x * 255.0f);  // R
        pixels[idx + 1] = (unsigned char)fminf(255.0f, finalColor.y * 255.0f);  // G
        pixels[idx + 2] = (unsigned char)fminf(255.0f, finalColor.z * 255.0f);  // B
        pixels[idx + 3] = 255;
    }
}

// ============== Renderer State ==============

struct CudaRenderer {
    // DX11 resources
    ID3D11Device* device;
    ID3D11DeviceContext* context;
    ID3D11Texture2D* texture;
    ID3D11Texture2D* stagingTexture;
    ID3D11ShaderResourceView* srv;

    // CUDA device memory
    unsigned char* d_pixels;
    unsigned char* h_pixels;  // Pinned host memory

    // Device memory for molecule data
    Atom* d_atoms;
    Bond* d_bonds;

    // Device memory for text
    char* d_text;

    // Dimensions
    int width;
    int height;
};

// ============== Renderer Implementation ==============

static bool createRenderTarget(CudaRenderer* r, int width, int height) {
    // Create DX11 texture for display
    D3D11_TEXTURE2D_DESC desc = {};
    desc.Width = width;
    desc.Height = height;
    desc.MipLevels = 1;
    desc.ArraySize = 1;
    desc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    desc.SampleDesc.Count = 1;
    desc.Usage = D3D11_USAGE_DEFAULT;
    desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;

    HRESULT hr = r->device->CreateTexture2D(&desc, nullptr, &r->texture);
    if (FAILED(hr)) {
        printf("Failed to create DX11 texture (hr=0x%08X)\n", hr);
        return false;
    }

    // Create staging texture for CPU->GPU copy
    desc.Usage = D3D11_USAGE_STAGING;
    desc.BindFlags = 0;
    desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

    hr = r->device->CreateTexture2D(&desc, nullptr, &r->stagingTexture);
    if (FAILED(hr)) {
        printf("Failed to create staging texture (hr=0x%08X)\n", hr);
        return false;
    }

    // Create shader resource view
    D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
    srvDesc.Format = desc.Format;
    srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
    srvDesc.Texture2D.MipLevels = 1;

    hr = r->device->CreateShaderResourceView(r->texture, &srvDesc, &r->srv);
    if (FAILED(hr)) {
        printf("Failed to create shader resource view (hr=0x%08X)\n", hr);
        return false;
    }

    // Allocate CUDA memory
    cudaMalloc(&r->d_pixels, width * height * 4);
    cudaMallocHost(&r->h_pixels, width * height * 4);

    r->width = width;
    r->height = height;

    printf("Render target created: %dx%d\n", width, height);
    return true;
}

static void destroyRenderTarget(CudaRenderer* r) {
    if (r->d_pixels) { cudaFree(r->d_pixels); r->d_pixels = nullptr; }
    if (r->h_pixels) { cudaFreeHost(r->h_pixels); r->h_pixels = nullptr; }
    if (r->srv) { r->srv->Release(); r->srv = nullptr; }
    if (r->stagingTexture) { r->stagingTexture->Release(); r->stagingTexture = nullptr; }
    if (r->texture) { r->texture->Release(); r->texture = nullptr; }
}

CudaRenderer* renderer_init(ID3D11Device* device, int width, int height) {
    // Get device context
    ID3D11DeviceContext* context;
    device->GetImmediateContext(&context);

    CudaRenderer* r = new CudaRenderer();
    memset(r, 0, sizeof(CudaRenderer));
    r->device = device;
    r->context = context;

    // Create render target
    if (!createRenderTarget(r, width, height)) {
        context->Release();
        delete r;
        return nullptr;
    }

    // Allocate device memory for molecule
    cudaMalloc(&r->d_atoms, MAX_ATOMS * sizeof(Atom));
    cudaMalloc(&r->d_bonds, MAX_BONDS * sizeof(Bond));
    cudaMalloc(&r->d_text, 256);  // Max text length

    printf("CUDA Renderer initialized (%dx%d)\n", width, height);
    return r;
}

void renderer_cleanup(CudaRenderer* r) {
    if (!r) return;

    destroyRenderTarget(r);

    if (r->d_atoms) cudaFree(r->d_atoms);
    if (r->d_bonds) cudaFree(r->d_bonds);
    if (r->d_text) cudaFree(r->d_text);
    if (r->context) r->context->Release();

    delete r;
    printf("CUDA Renderer cleaned up\n");
}

bool renderer_resize(CudaRenderer* r, int width, int height) {
    if (!r || (r->width == width && r->height == height)) return true;

    destroyRenderTarget(r);
    return createRenderTarget(r, width, height);
}

void renderer_render(CudaRenderer* r, const Molecule* mol, float rotX, float rotY, float zoom, float offsetX, float offsetY) {
    if (!r || !mol) return;

    // Copy molecule data to device
    cudaMemcpy(r->d_atoms, mol->atoms, mol->numAtoms * sizeof(Atom), cudaMemcpyHostToDevice);
    cudaMemcpy(r->d_bonds, mol->bonds, mol->numBonds * sizeof(Bond), cudaMemcpyHostToDevice);

    // Launch kernels
    dim3 blockSize(16, 16);
    dim3 gridSize((r->width + 15) / 16, (r->height + 15) / 16);

    clearKernel<<<gridSize, blockSize>>>(r->d_pixels, r->width, r->height);

    renderMoleculeKernel<<<gridSize, blockSize>>>(
        r->d_pixels, r->width, r->height,
        r->d_atoms, mol->numAtoms,
        r->d_bonds, mol->numBonds,
        rotX, rotY, zoom, offsetX, offsetY);

    cudaDeviceSynchronize();

    // Copy pixels from device to host
    cudaMemcpy(r->h_pixels, r->d_pixels, r->width * r->height * 4, cudaMemcpyDeviceToHost);

    // Copy to staging texture
    D3D11_MAPPED_SUBRESOURCE mapped;
    HRESULT hr = r->context->Map(r->stagingTexture, 0, D3D11_MAP_WRITE, 0, &mapped);
    if (SUCCEEDED(hr)) {
        // Copy row by row (handle pitch differences)
        for (int y = 0; y < r->height; y++) {
            memcpy((unsigned char*)mapped.pData + y * mapped.RowPitch,
                   r->h_pixels + y * r->width * 4,
                   r->width * 4);
        }
        r->context->Unmap(r->stagingTexture, 0);

        // Copy staging to final texture
        r->context->CopyResource(r->texture, r->stagingTexture);
    }
}

ID3D11ShaderResourceView* renderer_get_texture(CudaRenderer* r) {
    return r ? r->srv : nullptr;
}

void renderer_get_size(CudaRenderer* r, int* width, int* height) {
    if (r) {
        if (width) *width = r->width;
        if (height) *height = r->height;
    }
}

void renderer_render_text(CudaRenderer* r, const char* text, int x, int y, int scale) {
    if (!r || !text) return;

    int textLen = 0;
    while (textLen < 255 && text[textLen] != '\0') ++textLen;
    if (textLen == 0 || textLen > 255) return;

    // Copy text to device
    cudaMemcpy(r->d_text, text, textLen + 1, cudaMemcpyHostToDevice);

    // Launch text kernel
    dim3 blockSize(16, 16);
    dim3 gridSize((r->width + 15) / 16, (r->height + 15) / 16);

    renderTextKernel<<<gridSize, blockSize>>>(
        r->d_pixels, r->width, r->height,
        r->d_text, textLen, x, y, scale);

    cudaDeviceSynchronize();

    // Update the texture with new text
    cudaMemcpy(r->h_pixels, r->d_pixels, r->width * r->height * 4, cudaMemcpyDeviceToHost);

    D3D11_MAPPED_SUBRESOURCE mapped;
    HRESULT hr = r->context->Map(r->stagingTexture, 0, D3D11_MAP_WRITE, 0, &mapped);
    if (SUCCEEDED(hr)) {
        for (int row = 0; row < r->height; row++) {
            memcpy((unsigned char*)mapped.pData + row * mapped.RowPitch,
                   r->h_pixels + row * r->width * 4,
                   r->width * 4);
        }
        r->context->Unmap(r->stagingTexture, 0);
        r->context->CopyResource(r->texture, r->stagingTexture);
    }
}
