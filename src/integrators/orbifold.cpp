#include <mitsuba/core/orbifold.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/util.h>
#include <string>
#include <algorithm>


MTS_NAMESPACE_BEGIN

/*
 * Custom math functions which perform better than the ones provided by the standard library
 * max, floor, and ceil all perform expensive branching, these ones don't branch and give
 * a significant speedup to functions in this module which are called for every ray
 *
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
 *
 */
void OrbifoldData::InitX333OrbifoldData() {
    attenuate_callback = &OrbifoldData::attenuateX333;
    //    collapse_callback = &OrbifoldData::collapseX333;

    const double SQRT_3_OVER_2 = 0.8660254037844386;

    m_directionInfo[0].h_sep = m_scale.x;
    m_directionInfo[0].v_sep = SQRT_3_OVER_2 * m_scale.x;
    m_directionInfo[0].offset_array[0] = 0;
    m_directionInfo[0].offset_array[1] = 1.5 * m_scale.x;
    m_directionInfo[0].base = Vector2(0.0, 0.0);
}

void OrbifoldData::InitX2222OrbifoldData() {
    attenuate_callback = &OrbifoldData::attenuateX2222;
    m_directionInfo[0].v_sep = m_scale.z;
    m_directionInfo[1].v_sep = m_scale.x;
    m_directionInfo[0].base = Vector2(0.0, 0.0);
}

void OrbifoldData::InitXXOrbifoldData() {
    attenuate_callback = &OrbifoldData::attenuateXX;
    m_directionInfo[0].v_sep = m_scale.z;
    m_directionInfo[0].base = Vector2(0.0, 0.0);
}

void OrbifoldData::InitX632OrbifoldData() {
    attenuate_callback = &OrbifoldData::attenuateX632;
    //  collapse_callback = &OrbifoldData::collapseX632;

    const double SQRT_3_OVER_2 = 0.8660254037844386;

    m_directionInfo[0].h_sep = m_scale.x;
    m_directionInfo[0].v_sep = SQRT_3_OVER_2 * m_scale.x;
    m_directionInfo[0].offset_array[0] = 0;
    m_directionInfo[0].offset_array[1] = 1.5 * m_scale.x;
    m_directionInfo[0].base = Vector2(0.0, SQRT_3_OVER_2) * m_scale.x;

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
 * The result is stored in mirror_count or returned if there is only one mirror type intersected.
 *
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
 *
 */
Spectrum OrbifoldData::attenuateXX(const Ray& ray, const Intersection& isect) const {
    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const StaticDirectionInfo sdi(Vector2(0, 1), Vector2(1, 0), 0, 1, 0);

    const Vector3 nd(0, 0, 1);
    const Vector3 dn = normalize(ray.d);
    const Float cosThetaI = std::max(dot(nd, dn), dot(nd, -dn));

    Spectrum ret(1.0);

    unsigned mirrorCount[2] = {0, 0};
    countMirrorsHomogeneous(S, E, dir, sdi, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 2; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI, m_eta[i], m_k[i]);
        }
    }

    return ret;
}

Spectrum OrbifoldData::attenuateX2222(const Ray& ray, const Intersection& isect) const {
    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0), 0, 1, 0);
    const StaticDirectionInfo sdi2(Vector2(1, 0), Vector2(0, 1), 2, 3, 0);

    const Vector3 nd1(0, 0, 1), nd2(1, 0, 0);
    const Vector3 dn = normalize(ray.d);
    const Float cosThetaI1 = std::max(dot(nd1, dn), dot(nd1, -dn)), cosThetaI2 = std::max(dot(nd2, dn), dot(nd2, -dn));


    Spectrum ret(1.0);

    unsigned mirrorCount[4] = {0, 0, 0, 0};
    countMirrorsHomogeneous(S, E, dir, sdi1, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 4; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI1, m_eta[i], m_k[i]);
        }
    }

    memset(mirrorCount, 0, sizeof(mirrorCount));
    countMirrorsHomogeneous(S, E, dir, sdi2, m_directionInfo[1], mirrorCount);
    for (unsigned i = 0; i < 4; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI2, m_eta[i], m_k[i]);
        }
    }

    return ret;
}

Spectrum OrbifoldData::attenuateX442(const Ray& ray, const Intersection& isect) const {
    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const double ONE_OVER_SQRT_2 = 0.7071067811865475;

    const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0), 0, 1, 0);
    const StaticDirectionInfo sdi2(Vector2(1, 0), Vector2(0, 1), 1, 0, 0);
    const StaticDirectionInfo sdi3(Vector2(ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), Vector2(-ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), 2, 0, 0);
    const StaticDirectionInfo sdi4(Vector2(-ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), Vector2(ONE_OVER_SQRT_2, ONE_OVER_SQRT_2), 2, 0, 0);

    const Vector3 nd1(0, 0, 1),  nd2(1, 0, 0), nd3(ONE_OVER_SQRT_2, 0.0, ONE_OVER_SQRT_2), nd4(-ONE_OVER_SQRT_2, 0.0, ONE_OVER_SQRT_2);
    const Vector3 dn = normalize(ray.d);
    const Float cosThetaI1 = std::max(dot(nd1, dn), dot(nd1, -dn)),
                cosThetaI2 = std::max(dot(nd2, dn), dot(nd2, -dn)),
                cosThetaI3 = std::max(dot(nd3, dn), dot(nd3, -dn)),
                cosThetaI4 = std::max(dot(nd3, dn), dot(nd4, -dn));


    Spectrum ret(1.0);

    unsigned mirrorCount[3] = {0, 0, 0};
    countMirrorsHomogeneous(S, E, dir, sdi1, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI1, m_eta[i], m_k[i]);
        }
    }

    memset(mirrorCount, 0, sizeof(mirrorCount));
    countMirrorsHomogeneous(S, E, dir, sdi2, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI2, m_eta[i], m_k[i]);
        }
    }

    unsigned m = countMirrorsHomogeneous(S, E, dir, sdi3, m_directionInfo[1]);
    for (unsigned i = 0; i < m; i++) { ret *= fresnelConductorExact(cosThetaI3, m_eta[1], m_k[1]); }

    m = countMirrorsHomogeneous(S, E, dir, sdi4, m_directionInfo[1]);
    for (unsigned i = 0; i < m; i++) { ret *= fresnelConductorExact(cosThetaI4, m_eta[1], m_k[1]); }

    return ret;
}

Spectrum OrbifoldData::attenuateX333(const Ray& ray, const Intersection& isect) const {
    const double SQRT_3_OVER_2 = 0.8660254037844386;

    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0),  0, 1, 2);

    const StaticDirectionInfo sdi2(Vector2(-SQRT_3_OVER_2, 0.5), Vector2(0.5, SQRT_3_OVER_2), 2, 1, 0);

    const StaticDirectionInfo sdi3(Vector2(SQRT_3_OVER_2, 0.5), Vector2(-0.5, SQRT_3_OVER_2), 0, 1, 2);

    const Vector3 nd1(0, 0, 1), nd2(-SQRT_3_OVER_2, 0.0, 0.5), nd3(SQRT_3_OVER_2, 0.0, 0.5);
    const Vector3 dn = normalize(ray.d);
    const Float cosThetaI1 = std::max(dot(nd1, dn), dot(nd1, -dn)),
                cosThetaI2 = std::max(dot(nd2, dn), dot(nd2, -dn)),
                cosThetaI3 = std::max(dot(nd3, dn), dot(nd3, -dn));


    Spectrum ret(1.0);

    unsigned mirrorCount[3] = {0, 0, 0};
    countMirrorsHeterogeneous(S, E, dir, sdi1, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI1, m_eta[i], m_k[i]);
        }
    }

    memset(mirrorCount, 0, sizeof(mirrorCount));
    countMirrorsHeterogeneous(S, E, dir, sdi2, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI2, m_eta[i], m_k[i]);
        }
    }

    memset(mirrorCount, 0, sizeof(mirrorCount));
    countMirrorsHeterogeneous(S, E, dir, sdi3, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI3, m_eta[i], m_k[i]);
        }
    }

    return ret;
}

Spectrum OrbifoldData::attenuateX632(const Ray& ray, const Intersection& isect) const {
    const double SQRT_3_OVER_2 = 0.8660254037844386;

    const Vector2 S = Vector2(ray.o.x, ray.o.z) - m_directionInfo[0].base;
    const Vector2 E = Vector2(isect.p.x, isect.p.z) - m_directionInfo[0].base;
    const Vector2 dir = Vector2(ray.d.x, ray.d.z);

    const StaticDirectionInfo sdi1(Vector2(0, 1), Vector2(1, 0), 2, 0, 2);
    const StaticDirectionInfo sdi2(Vector2(-0.5, SQRT_3_OVER_2), Vector2(SQRT_3_OVER_2, 0.5), 1, 1, 1);
    const StaticDirectionInfo sdi3(Vector2(-SQRT_3_OVER_2, 0.5), Vector2(0.5, SQRT_3_OVER_2), 2, 0, 2);
    const StaticDirectionInfo sdi4(Vector2(-1, 0), Vector2(0, 1), 1, 1, 1);
    const StaticDirectionInfo sdi5(Vector2(SQRT_3_OVER_2, 0.5), Vector2(0.5, -SQRT_3_OVER_2), 2, 0, 2);
    const StaticDirectionInfo sdi6(Vector2(0.5, SQRT_3_OVER_2), Vector2(SQRT_3_OVER_2, -0.5), 1, 1, 1);


    const Vector3 nd1(0, 0, 1),  nd2(-0.5, 0.0, SQRT_3_OVER_2), nd3(-SQRT_3_OVER_2, 0.0, 0.5),
                  nd4(-1, 0, 0), nd5(SQRT_3_OVER_2, 0.0, 0.5),  nd6(0.5, 0.0, SQRT_3_OVER_2);
    const Vector3 dn = normalize(ray.d);
    const Float cosThetaI1 = std::max(dot(nd1, dn), dot(nd1, -dn)),
                cosThetaI2 = std::max(dot(nd2, dn), dot(nd2, -dn)),
                cosThetaI3 = std::max(dot(nd3, dn), dot(nd3, -dn)),
                cosThetaI4 = std::max(dot(nd4, dn), dot(nd4, -dn)),
                cosThetaI5 = std::max(dot(nd5, dn), dot(nd5, -dn)),
                cosThetaI6 = std::max(dot(nd6, dn), dot(nd6, -dn));

    Spectrum ret(1.0);

    unsigned m = countMirrorsHomogeneous(S, E, dir, sdi2, m_directionInfo[1]);
    for (unsigned i = 0; i < m; i++) { ret *= fresnelConductorExact(cosThetaI2, m_eta[1], m_k[1]); }

    m = countMirrorsHomogeneous(S, E, dir, sdi4, m_directionInfo[1]);
    for (unsigned i = 0; i < m; i++) { ret *= fresnelConductorExact(cosThetaI4, m_eta[1], m_k[1]); }

    m = countMirrorsHomogeneous(S, E, dir, sdi6, m_directionInfo[1]);
    for (unsigned i = 0; i < m; i++) { ret *= fresnelConductorExact(cosThetaI6, m_eta[1], m_k[1]); }

    unsigned mirrorCount[3] = {0, 0, 0};
    countMirrorsHeterogeneous(S, E, dir, sdi1, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI1, m_eta[i], m_k[i]);
        }
    }

    memset(mirrorCount, 0, sizeof(mirrorCount));
    countMirrorsHeterogeneous(S, E, dir, sdi3, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI3, m_eta[i], m_k[i]);
        }
    }

    memset(mirrorCount, 0, sizeof(mirrorCount));
    countMirrorsHeterogeneous(S, E, dir, sdi5, m_directionInfo[0], mirrorCount);
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < mirrorCount[i]; j++) {
            ret *= fresnelConductorExact(cosThetaI5, m_eta[i], m_k[i]);
        }
    }

    return ret;
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
