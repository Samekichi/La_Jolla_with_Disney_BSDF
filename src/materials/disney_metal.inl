#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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

    // Evaluate key parameters for DisneyMetal
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

	// Compute DisneyMetal's subcomponents
	// 0. common vars
	Vector3 h = normalize(dir_in + dir_out);
    Vector3 h_local = to_local(frame, h);
	Real n_dot_in = dot(frame.n, dir_in);  // cos(theta_in)
	Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
    // 1. F_m
	Spectrum F_m = schlick_fresnel(base_color, h_dot_out);
    //Spectrum F_m = base_color + (1 - base_color) * pow(1 - h_dot_out, 5);
    // 2. D_m
	Real D_m = GTR2(h_local, roughness, anisotropic);
    // 3. G_m
	Real G_m = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
               smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);
    // Compute the final result
    Spectrum f_metal = F_m * D_m * G_m / (4 * abs(n_dot_in));
    return f_metal;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    Vector3 h = normalize(dir_in + dir_out);
    Vector3 h_local = to_local(frame, h);
	Real n_dot_in = dot(frame.n, dir_in);

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // For metal, we use the ellipsoidal sampling from Heitz 2018
    // "Sampling the GGX Distribution of Visible Normals"
    // https://jcgt.org/published/0007/04/01/
    // this importance samples smith_masking(cos_theta_in) * GTR2(cos_theta_h, roughness) * cos_theta_out
    Real G_m = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic);
    Real D_m = GTR2(h_local, roughness, anisotropic);
    // (4 * cos_theta_v) is the Jacobian of the reflectiokn
    return (G_m * D_m) / (4 * n_dot_in);
    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    
    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real alpha_min = 0.0001;
    Real alpha_x = max(alpha_min, pow(roughness, 2) / aspect);
    Real alpha_y = max(alpha_min, pow(roughness, 2) * aspect);
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

    // Transform the micro normal to world space
    Vector3 h = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
