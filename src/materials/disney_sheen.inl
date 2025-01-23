#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    // Evaluate key parameters for DisneySheen
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	sheen_tint = std::clamp(sheen_tint, Real(0.01), Real(1));

    // Compute DisneySheen's subcomponents
    // 0. common vars
    Vector3 h = normalize(dir_in + dir_out);
    Real n_dot_out = dot(frame.n, dir_out);  // cos(theta_out)
    Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
    // 1. C_tint
    Spectrum C_tint = luminance(base_color) > 0 ? base_color / luminance(base_color) : make_const_spectrum(1);
    // 2. C_sheen
    Spectrum C_sheen = (1 - sheen_tint) + sheen_tint * C_tint;
    // 3. f_sheen
    Spectrum f_sheen = C_sheen * pow(1 - abs(h_dot_out), 5) * abs(n_dot_out);
    return f_sheen;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    // Same as Lambertian & DisneyDiffuse, we importance sample the cosine hemisphere domain.
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    // Same as Lambertian & DisneyDiffuse, we importance sample the cosine hemisphere domain.
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */ };
    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
