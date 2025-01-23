#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    
	// Evaluate key parameters for DisneyClearcoat
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

	// Compute DisneyClearcoat's subcomponents
	// 0. common vars
    Vector3 h = normalize(dir_in + dir_out);
	Vector3 h_local = to_local(frame, h);
	Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
	Real n_dot_in = dot(frame.n, dir_in);  // cos(theta_in)
    // 1. F_c
	Real R_0 = pow(1.5 - 1, 2) / pow(1.5 + 1, 2);  // R_0(eta), where we hardcoded eta = 1.5
	Real F_c = R_0 + (1 - R_0) * pow(1 - abs(h_dot_out), 5);
    // 2. D_c
    Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real D_c = ((alpha_g * alpha_g) - 1) / (c_PI * log(alpha_g * alpha_g) * (1 + (alpha_g * alpha_g - 1) * (h_local.z * h_local.z)));
    // 3. G_c
	/*// roughness = 0.25 -> ad-hoc fit, no clear geometric; this func pow(4) the 0.25, so need to sqrt
    Real G_c = smith_masking_gtr2(to_local(frame, dir_in), sqrt(0.25)) *
	        smith_masking_gtr2(to_local(frame, dir_out), sqrt(0.25)); */
    Real G_c = smith_masking_gtr2_alpha(to_local(frame, dir_in), 0.25) *
               smith_masking_gtr2_alpha(to_local(frame, dir_out), 0.25);  // roughness = 0.25 -> ad-hoc fit, no clear geometric
	// Compute the final result
	Spectrum f_clearcoat = make_const_spectrum( F_c * D_c * G_c / (4 * abs(n_dot_in)) );
    return f_clearcoat;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    
    // Evaluate key parameters for DisneyClearcoat
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Compute subcomponents
    // 0. common vars
    Vector3 h = normalize(dir_in + dir_out);
    Vector3 h_local = to_local(frame, h);
    Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
	Real n_dot_h = dot(frame.n, h);  // cos(theta_half)
    // 1. D_c
    Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real D_c = ((alpha_g * alpha_g) - 1) / (c_PI * log(alpha_g * alpha_g) * (1 + (alpha_g * alpha_g - 1) * (h_local.z * h_local.z)));
    // Compute final PDF
    return D_c * abs(n_dot_h) / (4 * abs(h_dot_out));
    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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

    // Evaluate key parameters for DisneyClearcoat
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
	// Compute subcomponents
    Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;

	// Compute parameters from uniform r.v.'s
    Real h_azimuth = c_TWOPI * rnd_param_uv.y;
    Real cos_h_elevation = sqrt((1 - pow(alpha_g * alpha_g, 1 - rnd_param_uv.x)) / (1 - alpha_g * alpha_g));
    Real sin_h_elevation = sqrt(1 - cos_h_elevation * cos_h_elevation);

    // Compute final h
    Vector3 h_local{sin_h_elevation * cos(h_azimuth), sin_h_elevation * sin(h_azimuth), cos_h_elevation};
	Vector3 h = to_world(frame, h_local);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, Real(1) /* roughness */ };
    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
