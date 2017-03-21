#include "orbifold.h"
#include <mitsuba/core/logger.h>
#include <string>
#include <algorithm>


MTS_NAMESPACE_BEGIN

/*
 * Custom math functions which perform better than the ones provided by the standard library
 * max, floor, and ceil all perform expensive branching, these ones don't branch and give
 * a significant speedup to functions in this module which are called for every ray
 */

static inline unsigned fastMax(unsigned x, unsigned y) {
  return (x ^ ((x ^ y) & -(x < y)));
}

static inline int fastFloor(Float x) {
  int i = (int)x; /* truncate */
  return i - ( i > x ); /* convert trunc to floor */
}

static inline int fastCeil(Float x) {
  int i = (int)x; /* truncate and add 1 */
  return i + 1 - ( i > x ); /* convert trunc to ceil */
}

template <typename T>
static inline int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

static inline void reflect(const Vector2& n, Vector3& v1, Vector3& v2) {
  const Float a00 = 1 - 2 * n[0] * n[0];
  const Float a01 = -2 * n[0] * n[1];
  const Float a10 = a01;
  const Float a11 = 1 - 2 * n[1] * n[1];

  const Float x1p = a00*v1.x + a01*v1.z;
  const Float z1p = a10*v1.x + a11*v1.z;

  const Float x2p = a00*v2.x + a01*v2.z;
  const Float z2p = a10*v2.x + a11*v2.z;

  v1.x = x1p;
  v1.z = z1p;
  v2.x = x2p;
  v2.z = z2p;
}

static inline double fastPow(const double base, const size_t exp) {
    double ret = 1.0;
    for (size_t i = 0; i < exp; i++) {
        ret *= base;
    }
    return ret;
}



bool OrbifoldData::setType(OrbifoldType type) {
    switch (type) {
    case OrbifoldType::OT_XX:
        InitXXOrbifoldData();
        break;
    case OrbifoldType::OT_X2222:
        InitX2222OrbifoldData();
        break;
    case OrbifoldType::OT_X442:
        InitX442OrbifoldData();
        break;
    case OrbifoldType::OT_X333:
        InitX333OrbifoldData();
        break;
    case OrbifoldType::OT_X632:
        InitX632OrbifoldData();
        break;
    case OrbifoldType::OT_NONE:
        break;
    case OrbifoldType::OT_ERROR:
        SLog(EError, "Orbifold type was not valid!");
    }

    return true;
}

/*
 * OrbifoldData initialization functions
 */

void OrbifoldData::InitX333OrbifoldData() {
    attenuate_callback = &OrbifoldData::attenuateX333;
//    collapse_callback = &OrbifoldData::collapseX333;

    const double SQRT_3_OVER_2 = 0.8660254037844386;
    const double SQRT_3_OVER_4 = 0.4330127018922193;
    const double SQRT_3_OVER_6 = 0.28867513459481287;

    m_directionInfo[0].h_sep = m_scale.x;
    m_directionInfo[0].v_sep = SQRT_3_OVER_2 * m_scale.x;
    m_directionInfo[0].offset_array[0] = 0;
    m_directionInfo[0].offset_array[1] = 1.5 * m_scale.x;
    //  direction_info[0].base = Vector2(-0.5, -SQRT_3_OVER_4) * scale.x();
    m_directionInfo[0].base = Vector2(0.0, 0.0); // Vector2(-0.5, -SQRT_3_OVER_6) * m_scale.x();
}

void OrbifoldData::InitX2222OrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateX2222;
  m_directionInfo[0].v_sep = m_scale.z;
  m_directionInfo[1].v_sep = m_scale.x;
  m_directionInfo[0].base = Vector2(0.0, 0.0); // Vector2(-1.0* m_scale.x, -1.0 * m_scale.z);
}

void OrbifoldData::InitXXOrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateXX;
  m_directionInfo[0].v_sep = m_scale.z;
  m_directionInfo[0].base = Vector2(0.0, 0.0); // Vector2(-1.0* m_scale.x, -1.0 * m_scale.z);
}

void OrbifoldData::InitX632OrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateX632;
//  collapse_callback = &OrbifoldData::collapseX632;

  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;
  const double SQRT_3_OVER_6 = 0.28867513459481287;

  m_directionInfo[0].h_sep = m_scale.x;
  m_directionInfo[0].v_sep = SQRT_3_OVER_2 * m_scale.x;
  m_directionInfo[0].offset_array[0] = 0;
  m_directionInfo[0].offset_array[1] = 1.5 * m_scale.x;
  m_directionInfo[0].base = Vector2(0.0, SQRT_3_OVER_2) * m_scale.x; // Vector2(0.5 - (1.0/3.0), SQRT_3_OVER_2 - SQRT_3_OVER_6) * scale.x;

  m_directionInfo[1].v_sep = 1.5 * m_scale.x;
}

void OrbifoldData::InitX442OrbifoldData() {
  const double SQRT_2 = 1.4142135623730951;

  attenuate_callback = &OrbifoldData::attenuateX442;
  m_directionInfo[0].v_sep = m_scale.x;
  m_directionInfo[1].v_sep = m_scale.x * SQRT_2;
  m_directionInfo[0].base = Vector2(0.0, 0.0); // Vector2(2.0/3.0, -1.0/3.0) * m_scale.x;
}





/*
 * Counts which mirrors are intersected by the straight path starting at S, and ending at E
 * where the mirrors are configured according to s_dir and d_dir.
 * The result is stored in mirror_count.
 */
void OrbifoldData::countMirrorsHeterogeneous(const Vector2& S, const Vector2& E,
    const Vector2& dir, const StaticDirectionInfo& s_dir, const DynamicDirectionInfo& d_dir,
    unsigned* mirror_count) const {

  // The projection of the start position onto the direction perpendicular to the mirror plane
  const Float So = dot(S, s_dir.o);

  // The projection of the path direction onto the direction perpendicular to the mirror plane
  const Float dir_dot_o = dot(dir, s_dir.o);

  // Compute the index of the first and last mirror planes intersected by the path
  int m0 = 0, mE = 0;
  if(dir_dot_o > 0) {
    m0 = fastCeil(So / d_dir.v_sep);
    mE = fastFloor(dot(E, s_dir.o) / d_dir.v_sep);
  } else if (dir_dot_o < 0) {
    m0 = fastFloor(So / d_dir.v_sep);
    mE = fastCeil(dot(E, s_dir.o) / d_dir.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return;
  }

  // The distance along the path of the first and last mirror intersections
  const Float k0 = (m0 * d_dir.v_sep - So) / dir_dot_o;
  const Float kE = (mE * d_dir.v_sep - So) / dir_dot_o;

  // Bail out if the path doesn't cross a mirror
  if(kE < k0) { return; }

  // The total number of mirrors intersected by the path
  const unsigned num_mirrors = abs(mE - m0);

  // dK is the distance along the path between mirrors
  const Float  dK = (kE - k0) / fastMax(num_mirrors, 1);

  // k is the distance along the path of each mirror plane intersected
  Float k = k0;
  // i keeps track of the type of mirror boundary
  unsigned i = abs(m0) % DynamicDirectionInfo::OFFSET_ARRAY_LEN;

  for(unsigned j = 0; j <= num_mirrors; j++) {
    // Where the ray intersects the mirror line
    const Vector2 P = S + k * dir;

    // Project the intersection onto the mirror plane, and figure out which type of mirror was intersected
    const int l = fastFloor((dot(P, s_dir.m) + d_dir.offset_array[i]) / d_dir.h_sep);
    mirror_count[s_dir.index_array[abs(l % static_cast<int>(StaticDirectionInfo::NUM_MIRRORS))]] += 1;

    k += dK;
    i = (i + 1) % DynamicDirectionInfo::OFFSET_ARRAY_LEN;
  }
}

void OrbifoldData::countMirrorsHomogeneous(const Vector2& S, const Vector2& E,
                        const Vector2& dir, const StaticDirectionInfo& s_dir,
                        const DynamicDirectionInfo& d_dir, unsigned* mirrors) const {
  // The projection of the start position onto the direction perpendicular to the mirror plane
  const Float So = dot(S, s_dir.o);

  // The projection of the path direction onto the direction perpendicular to the mirror plane
  const Float dir_dot_o = dot(dir, s_dir.o);

  // Compute the index of the first and last mirror planes intersected by the path
  int m0 = 0, mE = 0;
  if(dir_dot_o > 0) {
    m0 = fastCeil(So / d_dir.v_sep);
    mE = fastFloor(dot(E, s_dir.o) / d_dir.v_sep);
  } else if (dir_dot_o < 0) {
    // Note: Here the sign is inversed so that m0 > mE if the ray path intersects no mirrors
    m0 = -fastFloor(So / d_dir.v_sep);
    mE = -fastCeil(dot(E, s_dir.o) / d_dir.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return;
  }

  // If the ray intersects no mirrors, bail out
  if(m0 > mE) { return; }

  // The index in the mirror array of the first mirror intersected
  const unsigned m0_type = m0 & 1;

  // The total number of mirrors intersected by the path
  const unsigned num_mirrors = abs(mE - m0) + 1;

  // Used to add one in odd cases
  const unsigned mirror_parity = num_mirrors & 1;

  mirrors[s_dir.index_array[m0_type]] += (num_mirrors/2) + mirror_parity;
  mirrors[s_dir.index_array[(m0_type + 1) % 2]] += (num_mirrors/2);
}

unsigned OrbifoldData::countMirrorsHomogeneous(const Vector2& S, const Vector2& E,
                        const Vector2& dir, const StaticDirectionInfo& s_dir,
                        const DynamicDirectionInfo& d_dir) const {
  // The projection of the start position onto the direction perpendicular to the mirror plane
  const Float So = dot(S, s_dir.o);

  // The projection of the path direction onto the direction perpendicular to the mirror plane
  const Float dir_dot_o = dot(dir, s_dir.o);

  // Compute the index of the first and last mirror planes intersected by the path
  int m0 = 0, mE = 0;
  if(dir_dot_o > 0) {
    m0 = fastCeil(So / d_dir.v_sep);
    mE = fastFloor(dot(E, s_dir.o) / d_dir.v_sep);
  } else if (dir_dot_o < 0) {
    // Note: Here the sign is inversed so that m0 > mE if the ray path intersects no mirrors
    m0 = -fastFloor(So / d_dir.v_sep);
    mE = -fastCeil(dot(E, s_dir.o) / d_dir.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return 0;
  }

  // If the ray intersects no mirrors, bail out
  if(m0 > mE) { return 0; }

  // The number of mirrors intersected
  return abs(mE - m0) + 1;
}





/*
 * Ray attenuation functions. The attenuate callback pointer (attenuate_callback) is set to one of these methods
 */

Float OrbifoldData::attenuateXX(const Ray& ray, const Intersection& isect) const {
    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const StaticDirectionInfo sdi(Vector2(0, 1), Vector2(1, 0), 0, 1, 0);

    unsigned mirrorCount[2] = {0, 0};

    countMirrorsHomogeneous(S, E, dir, sdi, m_directionInfo[0], mirrorCount);

    return fastPow(m_r1, mirrorCount[0]) * fastPow(m_r2, mirrorCount[1]);
}

Float OrbifoldData::attenuateX2222(const Ray& ray, const Intersection& isect) const {
  const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
  const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
  const Vector2 dir = Vector2(ray.d.x, ray.d.z);

  const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0), 0, 1, 0);
  const StaticDirectionInfo sdi2(Vector2(1, 0), Vector2(0, 1), 2, 3, 0);

  unsigned mirrorCount[4] = {0, 0, 0, 0};

  countMirrorsHomogeneous(S, E, dir, sdi1, m_directionInfo[0], mirrorCount);
  countMirrorsHomogeneous(S, E, dir, sdi2, m_directionInfo[1], mirrorCount);

  return fastPow(m_r1, mirrorCount[0]) * fastPow(m_r2, mirrorCount[1]) * fastPow(m_r3, mirrorCount[2]) * fastPow(m_r4, mirrorCount[3]);
}

Float OrbifoldData::attenuateX442(const Ray& ray, const Intersection& isect) const {
    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const double ONE_OVER_SQRT_2 = 0.7071067811865475;

    const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0), 0, 1, 0);
    const StaticDirectionInfo sdi2(Vector2(1, 0), Vector2(0, 1), 1, 0, 0);
    const StaticDirectionInfo sdi3(Vector2(ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), Vector2(-ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), 2, 0, 0);
    const StaticDirectionInfo sdi4(Vector2(-ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), Vector2(ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), 2, 0, 0);

    unsigned mirror_count[3] = {0, 0, 0};

    countMirrorsHomogeneous(S, E, dir, sdi1, m_directionInfo[0], mirror_count);
    countMirrorsHomogeneous(S, E, dir, sdi2, m_directionInfo[0], mirror_count);
    mirror_count[2] += countMirrorsHomogeneous(S, E, dir, sdi3, m_directionInfo[1]);
    mirror_count[2] += countMirrorsHomogeneous(S, E, dir, sdi4, m_directionInfo[1]);

    return fastPow(m_r1, mirror_count[0]) * fastPow(m_r2, mirror_count[1]) * fastPow(m_r3, mirror_count[2]);
}

Float OrbifoldData::attenuateX333(const Ray& ray, const Intersection& isect) const {
    const double SQRT_3_OVER_2 = 0.8660254037844386;
    const double SQRT_3_OVER_4 = 0.4330127018922193;

    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0),  0, 1, 2);

    const StaticDirectionInfo sdi2(Vector2(-SQRT_3_OVER_2, 0.5), Vector2(0.5, SQRT_3_OVER_2), 2, 1, 0);

    const StaticDirectionInfo sdi3(Vector2(SQRT_3_OVER_2, 0.5), Vector2(-0.5, SQRT_3_OVER_2), 0, 1, 2);

    unsigned mirrorCount[3] = {0, 0, 0};

    countMirrorsHeterogeneous(S, E, dir, sdi1, m_directionInfo[0], mirrorCount);
    countMirrorsHeterogeneous(S, E, dir, sdi2, m_directionInfo[0], mirrorCount);
    countMirrorsHeterogeneous(S, E, dir, sdi3, m_directionInfo[0], mirrorCount);

    return fastPow(m_r1, mirrorCount[0]) * fastPow(m_r2, mirrorCount[1]) * fastPow(m_r3, mirrorCount[2]);
}

Float OrbifoldData::attenuateX632(const Ray& ray, const Intersection& isect) const {
    const double SQRT_3_OVER_2 = 0.8660254037844386;
    const double SQRT_3_OVER_4 = 0.4330127018922193;

    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0), 2, 0, 2);

    // Homogeneous
    const StaticDirectionInfo sdi2(Vector2(-0.5, SQRT_3_OVER_2), Vector2(SQRT_3_OVER_2, 0.5), 1, 1, 1);

    const StaticDirectionInfo sdi3(Vector2(-SQRT_3_OVER_2, 0.5), Vector2(0.5, SQRT_3_OVER_2), 2, 0, 2);

    // Homogeneous
    const StaticDirectionInfo sdi4(Vector2(-1, 0), Vector2(0, 1), 1, 1, 1);


    const StaticDirectionInfo sdi5(Vector2(SQRT_3_OVER_2, 0.5), Vector2(0.5, -SQRT_3_OVER_2), 2, 0, 2);

    // Homogeneous
    const StaticDirectionInfo sdi6(Vector2(0.5, SQRT_3_OVER_2), Vector2(SQRT_3_OVER_2, -0.5), 1, 1, 1);

    unsigned mirror_count[3] = {0, 0, 0};

    mirror_count[1] = countMirrorsHomogeneous(S, E, dir, sdi2, m_directionInfo[1]) +
                      countMirrorsHomogeneous(S, E, dir, sdi4, m_directionInfo[1]) +
                      countMirrorsHomogeneous(S, E, dir, sdi6, m_directionInfo[1]);

    countMirrorsHeterogeneous(S, E, dir, sdi1, m_directionInfo[0], mirror_count);
    countMirrorsHeterogeneous(S, E, dir, sdi3, m_directionInfo[0], mirror_count);
    countMirrorsHeterogeneous(S, E, dir, sdi5, m_directionInfo[0], mirror_count);

    return fastPow(m_r1, mirror_count[0]) * fastPow(m_r2, mirror_count[1]) * fastPow(m_r3, mirror_count[2]);
}


OrbifoldData::OrbifoldType OrbifoldData::OrbifoldTypeFromString(const std::string& str) {
    std::string lstr = str;
    std::transform(lstr.begin(), lstr.end(), lstr.begin(), ::tolower);

    if (lstr == "**") {
        return OrbifoldType::OT_XX;
    } else if (lstr == "*2222") {
        return OrbifoldType::OT_X2222;
    } if (lstr == "*333") {
        return OrbifoldType::OT_X333;
    } if (lstr == "*442") {
        return OrbifoldType::OT_X442;
    } if (lstr == "*632") {
        return OrbifoldType::OT_X632;
    } else if (lstr == "none") {
        return OrbifoldType::OT_NONE;
    } else {
        return OrbifoldType::OT_ERROR;
    }

    return OrbifoldType::OT_ERROR;
}

std::string OrbifoldData::OrbifoldTypeToString(const OrbifoldData::OrbifoldType& ot) {
    switch(ot) {
    case OrbifoldType::OT_XX:
        return std::string("**");
    case OrbifoldType::OT_X2222:
        return std::string("*2222");
    case OrbifoldType::OT_X333:
        return std::string("*333");
    case OrbifoldType::OT_X442:
        return std::string("*442");
    case OrbifoldType::OT_X632:
        return std::string("*632");
    case OrbifoldType::OT_NONE:
        return std::string("none");
    case OrbifoldType::OT_ERROR:
        return std::string("error");
    }

    return std::string("error");
}


MTS_NAMESPACE_END
