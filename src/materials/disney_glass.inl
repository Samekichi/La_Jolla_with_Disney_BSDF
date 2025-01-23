#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw

	// Evaluate Key Parameters for DisneyGlass
	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

	// Compute DisneyGlass's subcomponents
	// 0. half-vector
    Vector3 h;
    if (reflect) {
		h = normalize(dir_in + dir_out);
	}
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }
    Real h_dot_in = dot(h, dir_in);  // cos(theta_half_in)
	Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
	Real n_dot_in = dot(frame.n, dir_in);  // cos(theta_in)
	Real n_dot_out = dot(frame.n, dir_out);  // cos(theta_out)
 //   // 1. R_s
 //   Real R_s = (h_dot_in - eta * h_dot_out) / (h_dot_in + eta * h_dot_out);
 //   // 2. R_p
	//Real R_p = (eta * h_dot_in - h_dot_out) / (eta * h_dot_in + h_dot_out);
 //   // 3. F_g
 //   Real F_g = (R_s * R_s + R_p * R_p) / 2;
    Real F_g = fresnel_dielectric(h_dot_in, eta);
    // 4. D_g
    Real D_g = GTR2(to_local(frame, h), roughness, anisotropic);
    // 5. G_g
    Real G_g = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
               smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);
    
    // Compute the final result
    Spectrum f_glass;
    if (reflect) {
        // Case 1: Reflection
		f_glass = base_color * F_g * D_g * G_g / (4 * abs(n_dot_in));
    }
    else {
        // Case 2: Refraction
        f_glass = sqrt(base_color) * (1 - F_g) * D_g * G_g * abs(h_dot_in * h_dot_out) / 
                  (abs(n_dot_in) * (h_dot_in + eta * h_dot_out) * (h_dot_in + eta * h_dot_out));
    }
    return f_glass;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw
    
    // Evaluate Key Parameters for DisneyGlass
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Clamp to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
	anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));

    // half-vector
    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    // Common subcomponents
	Real h_dot_in = dot(h, dir_in);  // cos(theta_half_in)
	Real n_dot_in = dot(frame.n, dir_in);  // cos(theta_in)
    // This importance samples smith_masking(cos_theta_in) * GTR2(cos_theta_h, roughness) * cos_theta_out;
    // Note that we use the incoming direction
    // for evaluating the Fresnel reflection amount.
    // We can also use outgoing direction -- then we would need to
    // use 1/bsdf.eta and we will get the same result.
    Real F_g = fresnel_dielectric(h_dot_in, eta);
    Real D_g = GTR2(to_local(frame, h), roughness, anisotropic);
    Real G_g = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
               smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);
    
    if (reflect) {
		return F_g * D_g * G_g / (4 * abs(n_dot_in));
    }
    else {
        Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
        return (1 - F_g) * D_g * G_g * abs(h_dot_in * h_dot_out) /
            (abs(n_dot_in) * (h_dot_in + eta * h_dot_out) * (h_dot_in + eta * h_dot_out));
    }

    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Evaluate Key Parameters for DisneyGlass
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Clamp to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real alpha_min = 0.0001;
    Real alpha_x = max(alpha_min, (roughness * roughness) / aspect);
    Real alpha_y = max(alpha_min, (roughness * roughness) * aspect);
    Vector3 local_micro_normal =
        sample_visible_normals(to_local(frame, dir_in), alpha_x, alpha_y, rnd_param_uv);

    // Transform the micro normal to world space
    Vector3 h = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(h, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        // Reflect over the world space normal
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, roughness /* roughness */
        };
    }
    else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            h = -h;
        }
        Real h_dot_out = sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * h;
        return BSDFSampleRecord{ refracted, eta, roughness };
    }
    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
