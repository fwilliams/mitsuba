#include <mitsuba/render/scene.h>


MTS_NAMESPACE_BEGIN

class DepthIntegrator : public SamplingIntegrator {
public:
    MTS_DECLARE_CLASS()

    DepthIntegrator(const Properties& props) : SamplingIntegrator(props) {
        Spectrum defaultColor;
        defaultColor.fromLinearRGB(1.0f, 1.0f, 1.0);
        m_color = defaultColor;
//        m_color = props.getSpectrum("color", defaultColor);
    }

    /// Unserialize from a binary data stream
    DepthIntegrator(Stream* stream, InstanceManager* instMgr) : SamplingIntegrator(stream, instMgr) {
        SamplingIntegrator::serialize(stream, instMgr);
        m_color = Spectrum(stream);
        m_maxDist = stream->readFloat();
    }

    /// Serialize to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const {
    SamplingIntegrator::serialize(stream, manager);
        m_color.serialize(stream);
        stream->writeFloat(m_maxDist);
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        if (rRec.rayIntersect(ray)) {
            Float distance = rRec.its.t;

            return Spectrum(1.0 - distance/m_maxDist) * m_color;
        }

        return Spectrum(0.0);
    }

    /// Preprocess function -- called on the initiating machine
    bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job, int sceneResID, int cameraResID, int samplerResID) {
        SamplingIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);

        const AABB& sceneAABB = scene->getAABB();

        Point cameraPos = scene->getSensor()->getWorldTransform()->eval(0).transformAffine(Point(0.0));
        for (int i = 0; i < 8; i++) {
            m_maxDist = std::max(m_maxDist, (cameraPos - sceneAABB.getCorner(i)).length());
        }

        return true;
    }


private:
    Spectrum m_color;
    Float m_maxDist;
};

MTS_IMPLEMENT_CLASS_S(DepthIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(DepthIntegrator, "Generate a depth image")

MTS_NAMESPACE_END
