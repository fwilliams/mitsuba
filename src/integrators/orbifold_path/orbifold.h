#pragma once
#ifndef __MITSUBA_CORE_ORBIFOLD_H
#define __MITSUBA_CORE_ORBIFOLD_H

#include <mitsuba/render/shape.h>


MTS_NAMESPACE_BEGIN

struct OrbifoldData {

    enum class OrbifoldType {
        OT_XX,
        OT_X2222,
        OT_X333,
        OT_X442,
        OT_X632,
        OT_NONE,
        OT_ERROR
    };

    OrbifoldData() :
        m_incompletenessMode(false),
        m_scale(Vector3(1, 1, 1)),
        m_numKernelTiles(1),
        m_type(OrbifoldType::OT_NONE),
        attenuate_callback(&OrbifoldData::attenuateTrivial),
        collapse_callback(&OrbifoldData::collapseTrivial) {}

    ///
    /// \brief m_incompletenessMode
    /// Whether or not we are rendering an incompleteness image
    bool m_incompletenessMode;

    ///
    /// A scaling factor for the scense
    ///
    Vector3 m_scale;

    ///
    /// The index of refraction and reflectivity of each mirror.
    /// Depending on the orbifold type, not all of these are used
    ///
    Spectrum m_k[4];
    Spectrum m_eta[4];

    ///
    /// The number of tiles in the kernel being rendered
    ///
    unsigned m_numKernelTiles;

    ///
    /// Sets the type of orbifold, which will affect which attenuate and collapse functions are used
    ///
    bool setType(OrbifoldType type);

    ///
    /// Computes the attenuation along a ray which may intersect one or more mirror boundaries
    /// in the unfolded scene
    ///
    inline Spectrum attenuate(const Ray& r, const Intersection& i) const {
        return (this->*attenuate_callback)(r, i);
    }

    ///
    /// Transforms (collapses) a point in the covering space into the fundamental domain.
    /// This function returns the transformed point. out_fixed_dir is mutated to point
    /// in the direction corresponding to the transformation
    ///
    inline Vector3 collapse(const Vector3& point, Vector3& out_fixed_dir) const {
        return (this->*collapse_callback)(point, out_fixed_dir);
    }

    ///
    /// Transforms (collapses) a point in the covering space into the fundamental domain.
    /// This function returns the transformed point.
    ///
    inline Vector3 collapse(const Vector3& point) const {
        Vector3 dummy;
        return (this->*collapse_callback)(point, dummy);
    }

    ///
    /// \brief OrbifoldTypeFromString
    /// \param str
    /// \return The orbifold type represented by the string
    ///
    static OrbifoldType OrbifoldTypeFromString(const std::string& str);

    ///
    /// \brief OrbifoldTypeToString
    /// \param ot
    /// \return The string representation of an orbifold type
    ///
    static std::string OrbifoldTypeToString(const OrbifoldType& ot);


private:

    ///
    /// Convert a kernel radius to a number of tiles for each group
    ///
    static inline unsigned X333_KernelRadToTileNum(unsigned kernel_rad) {
        unsigned ret = 1;
        for(unsigned i = 1; i <= kernel_rad; i++) {
            ret += 6 * i;
        }
        return ret * 6;
    }
    static inline unsigned X2222_KernelRadToTileNum(unsigned kernel_rad) {
        unsigned ret = 1;
        for(unsigned i = 1; i <= kernel_rad; i++) {
            ret += 8 * i;
        }
        return ret * 4;
    }
    static inline unsigned XX_KernelRadToTileNum(unsigned kernel_rad) {
        return 4 * kernel_rad + 1;
    }
    static inline unsigned X442_KernelRadToTileNum(unsigned kernel_rad) {
        return X2222_KernelRadToTileNum(kernel_rad) * 2;
    }
    static inline unsigned X632_KernelRadToTileNum(unsigned kernel_rad) {
        return X333_KernelRadToTileNum(kernel_rad) * 2;
    }

    ///
    /// Static information about mirror directions
    ///
    struct StaticDirectionInfo {
        StaticDirectionInfo(const Vector2& oo, const Vector2& mm, unsigned ia1, unsigned ia2, unsigned ia3) :
            o(oo), m(mm), index_array{ia1, ia2, ia3} {}
        Vector2 o, m;
        unsigned index_array[3];
        static const unsigned NUM_MIRRORS = 3;
    };

    ///
    /// Runtime information about mirror directions
    ///
    struct DynamicDirectionInfo {
        Vector2 base;
        Float v_sep, h_sep;
        Float offset_array[2];
        static const unsigned OFFSET_ARRAY_LEN = 2;
    };

    ///
    /// Mirror attenuation callback type
    ///
    typedef Spectrum (OrbifoldData::*OrbifoldAttenuationCallback)(const Ray& r, const Intersection& isect) const;

    ///
    /// Point collapse callback (see collapse)
    ///
    typedef Vector3 (OrbifoldData::*OrbifoldCollapseCallback)(const Vector3& pt, Vector3& ofd) const;


    ///
    /// The type of orbifold to use
    ///
    OrbifoldType m_type;

    ///
    /// Information about adjacent mirror directions set at runtime
    ///
    DynamicDirectionInfo m_directionInfo[2];

    ///
    /// The callback function which calculates attenuation from mirror boundary intersections
    /// This gets set at runtime to one of the attenuate*() functions based on the type of orbifold
    /// specified.
    ///
    OrbifoldAttenuationCallback attenuate_callback;

    ///
    /// The callback function which collapses a virtual point into the fundamental domain.
    /// This gets set at runtime to one of the collapse*() functions based on the type of orbifold specified.
    ///
    OrbifoldCollapseCallback collapse_callback;

    ///
    /// Called when each mirror plane consists of a repeating pattern of different mirrors.
    /// Counts which mirrors are intersected by the straight path starting at S, and ending at E
    /// where the mirrors are configured according to s_dir and d_dir.
    /// The result is stored in mirror_count.
    ///
    void countMirrorsHeterogeneous(const Vector2& S, const Vector2& E,
                                   const Vector2& dir, const StaticDirectionInfo& s_dir,
                                   const DynamicDirectionInfo& d_dir, unsigned* mirrors) const;


    void countMirrorsHomogeneous(const Vector2& S, const Vector2& E,
                                 const Vector2& dir, const StaticDirectionInfo& s_dir,
                                 const DynamicDirectionInfo& d_dir, unsigned* mirrors) const;

    unsigned countMirrorsHomogeneous(const Vector2& S, const Vector2& E,
                                     const Vector2& dir, const StaticDirectionInfo& s_dir,
                                     const DynamicDirectionInfo& d_dir) const;

    Spectrum attenuateTrivial(const Ray& ray, const Intersection& isect) const {
        return Spectrum(1.0);
    }
    //  Float attenuateX333(const Ray& ray, const Intersection& isect) const;
    Spectrum attenuateX2222(const Ray& ray, const Intersection& isect) const;
    Spectrum attenuateXX(const Ray& ray, const Intersection& isect) const;
    Spectrum attenuateX333(const Ray& ray, const Intersection& isect) const;
    Spectrum attenuateX632(const Ray& ray, const Intersection& isect) const;
    Spectrum attenuateX442(const Ray& ray, const Intersection& isect) const;

    Vector3 collapseTrivial(const Vector3& pt, Vector3& o) const { return pt; }
    //  Vector3 collapseX333(const Vector3& pt, Vector3& o) const;
    //  Vector3 collapseX632(const Vector3& pt, Vector3& o) const;
    //  Vector3 collapseX442(const Vector3& pt, Vector3& o) const;
    //  Vector3 collapseX2222(const Vector3& pt, Vector3& o) const;
    //  Vector3 collapseXX(const Vector3& pt, Vector3& o) const;

    ///
    /// Initialize this structure based on the type of orbifold present
    ///
    void InitX333OrbifoldData();
    void InitX2222OrbifoldData();
    void InitXXOrbifoldData();
    void InitX442OrbifoldData();
    void InitX632OrbifoldData();
};


MTS_NAMESPACE_END

#endif // __MITSUBA_CORE_ORBIFOLD_H

