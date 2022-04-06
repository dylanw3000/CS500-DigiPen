
#include "geom.h"
#include "raytrace.h"
#include "acceleration.h"

#include <bvh/sweep_sah_builder.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/primitive_intersectors.hpp>

/////////////////////////////
// Vector and ray conversions
Ray RayFromBvh(const bvh::Ray<float> &r)
{
    return Ray(vec3FromBvh(r.origin), vec3FromBvh(r.direction));
}
bvh::Ray<float> RayToBvh(const Ray &r)
{
    return bvh::Ray<float>(vec3ToBvh(r.origin), vec3ToBvh(r.dir));
    // return bvh::Ray<float>(vec3ToBvh(r.o), vec3ToBvh(r.d));
}


/////////////////////////////
// SimpleBox
bvh::Vector3<float> vec3ToBvh(const vec3& v)
{
    return bvh::Vector3<float>(v[0],v[1],v[2]);
}

vec3 vec3FromBvh(const bvh::Vector3<float>& v)
{
    return vec3(v[0],v[1],v[2]);
}

SimpleBox::SimpleBox(): bvh::BoundingBox<float>() {}
SimpleBox::SimpleBox(const vec3 v): bvh::BoundingBox<float>(vec3ToBvh(v)) {}

SimpleBox& SimpleBox::extend(const vec3 v)
{
    bvh::BoundingBox<float>::extend(vec3ToBvh(v));
    return *this;
}


/////////////////////////////
// BvhShape

SimpleBox BvhShape::bounding_box() const
{
    //  Return the shape's bounding box.
    // return SimpleBox(); // FIX THIS
    // if(shape->accelerationBox != nullptr) return *shape->accelerationBox;
    // return SimpleBox(); // FIX THIS

    SimpleBox out;

    if (dynamic_cast<Sphere*>(shape) != nullptr) {
        Sphere* s = static_cast<Sphere*>(shape);
        
        out = SimpleBox(s->pos - vec3(1.f)*s->radius);
        out.extend(s->pos + vec3(1.f)*s->radius);
        
    }
    else if (dynamic_cast<Box*>(shape) != nullptr) {
        Box* b = static_cast<Box*>(shape);
        out = SimpleBox(b->b);
        out.extend(b->b + b->d);
    }
    else if (dynamic_cast<Cylinder*>(shape) != nullptr) {
        Cylinder* c = static_cast<Cylinder*>(shape);
        if (true) {
            out = SimpleBox(c->base);
            out.extend(c->base + vec3(c->radius));
            out.extend(c->base - vec3(c->radius));
            out.extend(c->base + c->ang);
            out.extend(c->base + c->ang + vec3(c->radius));
            out.extend(c->base + c->ang - vec3(c->radius));
        }
    }
    else if (dynamic_cast<Triangle*>(shape) != nullptr) {
        Triangle* t = static_cast<Triangle*>(shape);
        
        out = SimpleBox(t->v0);
        out.extend(t->v1);
        out.extend(t->v2);
        
    }

    return out;
}

bvh::Vector3<float> BvhShape::center() const
{
    return bounding_box().center();
}
    
std::optional<Intersection> BvhShape::intersect(const bvh::Ray<float>& bvhray) const
{
    // Intersect RayFromBvh(bvhray) with shape;  store result in an Intersection
    // If no intersection,
    //    return std::nullopt;
    // If intersection's t value < bvhray.tmin  or > bvhray.tmax
    //    return std::nullopt;
    // else return
    //    return the Intersection
    
    // return Intersection();  // FIX THIS 

    Intersection out = shape->Intersect(RayFromBvh(bvhray));
    if (out.collision) {
        if(out.t < bvhray.tmin || out.t > bvhray.tmax) return std::nullopt;
        return out;
    }
    return std::nullopt;
}

AccelerationBvh::AccelerationBvh(std::vector<Shape*> &objs)
{
    // Wrap all Shape*'s with a bvh specific instance
    for (Shape* shape:  objs) {
        shapeVector.emplace_back(shape); }

    // Magic found in the bvh examples:
    auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(shapeVector.data(),
                                                                     shapeVector.size());
    auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), shapeVector.size());

    bvh::SweepSahBuilder<bvh::Bvh<float>> builder(bvh);
    builder.build(global_bbox, bboxes.get(), centers.get(), shapeVector.size());
}

Intersection AccelerationBvh::intersect(const Ray& ray)
{
    bvh::Ray<float> bvhRay = RayToBvh(ray);

    // Magic found in the bvh examples:
    bvh::ClosestPrimitiveIntersector<bvh::Bvh<float>, BvhShape> intersector(bvh, shapeVector.data());
    bvh::SingleRayTraverser<bvh::Bvh<float>> traverser(bvh);

    auto hit = traverser.traverse(bvhRay, intersector);
    if (hit) {
        return hit->intersection;
    }
    else
        return  Intersection();  // Return an IntersectionRecord which indicates NO-INTERSECTION


}
