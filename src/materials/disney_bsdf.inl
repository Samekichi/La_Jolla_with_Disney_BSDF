#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    bool is_inside = dot(vertex.geometric_normal, dir_in) <= 0;

    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw

    // Evaluate key parameters for DisneyGlass
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);            
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    roughness = std::clamp(roughness, Real(0.01), Real(1));
    specular_transmission = std::clamp(specular_transmission, Real(0.01), Real(1));
    metallic = std::clamp(metallic, Real(0.01), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    subsurface = std::clamp(subsurface, Real(0.01), Real(1));
    specular = std::clamp(specular, Real(0.01), Real(1));
    specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
    sheen = std::clamp(sheen, Real(0.01), Real(1));
    sheen_tint = std::clamp(sheen_tint, Real(0.01), Real(1));
    clearcoat = std::clamp(clearcoat, Real(0.01), Real(1));
    clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));

    // Compute weight for each component
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_sheen = (1 - metallic) * sheen;
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;

    // Check if incident inside (inside = only glass is needed, other weights all 0)
    Spectrum f_glass = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    if (is_inside) {
        return w_glass * f_glass;
    }

    // Compute remaining DisneyBSDF's components
    Spectrum f_diffuse = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
	Spectrum f_clearcoat = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
	Spectrum f_sheen = operator()(DisneySheen{ bsdf.base_color, bsdf.sheen_tint });
	Spectrum f_metal;  // need modified Fresnel term `F_m` to include an achromatic specular component, with control para `specular_tint` to potentially make it closer to base_color
    // compute DisneyMetal
    // 0. common vars
    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
        //h = normalize(dir_in + dir_out);
    }
    // flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }
    Vector3 h_local = to_local(frame, h);
    Real n_dot_in = dot(frame.n, dir_in);  // cos(theta_in)
    Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
    // 1. F_m
    Spectrum C_tint = luminance(base_color) > 0 ? base_color / luminance(base_color) : make_const_spectrum(1);
    Spectrum K_s = (1 - specular_tint) + specular_tint * C_tint;
	Real R_0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    Spectrum C_0 = specular * R_0 * (1 - metallic) * K_s + metallic * base_color;
	Spectrum F_m = C_0 + (1 - C_0) * pow(1 - abs(h_dot_out), 5);
    // 2. D_m
    Real D_m = GTR2(h_local, roughness, anisotropic);
    // 3. G_m
    Real G_m = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
               smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);
    // final f_metal
    f_metal = F_m * D_m * G_m / (4 * abs(n_dot_in));

	// Compute the final result
    return w_diffuse * f_diffuse + 
		   w_sheen * f_sheen + 
           w_metal * f_metal + 
           w_clearcoat * f_clearcoat + 
           w_glass * f_glass;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    bool is_inside = dot(vertex.geometric_normal, dir_in) <= 0;
    
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw
    
    // Evaluate key parameters for DisneyGlass
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    specular_transmission = std::clamp(specular_transmission, Real(0.01), Real(1));
    metallic = std::clamp(metallic, Real(0.01), Real(1));
    subsurface = std::clamp(subsurface, Real(0.01), Real(1));
    specular = std::clamp(specular, Real(0.01), Real(1));
    specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
    clearcoat = std::clamp(clearcoat, Real(0.01), Real(1));
    clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));

    // Check if incident inside (inside = only glass is needed, other weights all 0)
    Real p_glass = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    if (is_inside) {
        return p_glass;
    }

    // Compute weight for each component
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;

    // Compute remaining DisneyBSDF's components (DisneySheen is negelected in pdf)
    Real p_diffuse = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    Real p_clearcoat = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
	//Real p_clearcoat = 0;
    // compute subcomponents
    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
        //h = normalize(dir_in + dir_out);
    }
    // flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }
    Vector3 h_local = to_local(frame, h);
    Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
    Real n_dot_h = dot(frame.n, h);  // cos(theta_half)
    // 1. D_c
    Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real D_c = ((alpha_g * alpha_g) - 1) / (c_PI * log(alpha_g * alpha_g) * (1 + (alpha_g * alpha_g - 1) * (h_local.z * h_local.z)));
    // Compute final PDF
    //p_clearcoat = D_c * abs(n_dot_h) / (4 * abs(h_dot_out));

    Real p_metal = 0;  // need to modify p_metal to consider different h
    Real n_dot_in = dot(frame.n, dir_in);  // cos(theta_in)
    // 1. D_m
    Real D_m = GTR2(h_local, roughness, anisotropic);
    // 2. G_m
    Real G_m = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
               smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);
    // final p_metal
    p_metal = D_m * G_m / (4 * abs(n_dot_in));
    
    // Normalize weights
    Real w = w_diffuse + w_metal + w_clearcoat + w_glass;
    w_diffuse /= w;
    w_metal /= w;
    w_clearcoat /= w;
    w_glass /= w;

	// Compute the final result
    return w_diffuse * p_diffuse +
           w_metal * p_metal +
           w_clearcoat * p_clearcoat +
           w_glass * p_glass;
    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool is_inside = dot(vertex.geometric_normal, dir_in) <= 0;

    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // Evaluate key parameters for DisneyBSDF
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    specular_transmission = std::clamp(specular_transmission, Real(0.01), Real(1));
    metallic = std::clamp(metallic, Real(0.01), Real(1));
    subsurface = std::clamp(subsurface, Real(0.01), Real(1));
    specular = std::clamp(specular, Real(0.01), Real(1));
    specular_tint = std::clamp(specular_tint, Real(0.01), Real(1));
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    clearcoat = std::clamp(clearcoat, Real(0.01), Real(1));
    clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));

	// Check if incident inside (inside = only glass is needed, other weights all 0)
    if (is_inside) {
        return operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    }

    // Compute raw weight for each component
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
    
    // Normalize weights
    Real w = w_diffuse + w_metal + w_clearcoat + w_glass;
    w_diffuse /= w;
	w_metal /= w;
	w_clearcoat /= w;
	w_glass /= w;

    // Importance sampling
	Real rnd = rnd_param_w;
    std::optional<BSDFSampleRecord> sample;
 //   if (rnd < w_diffuse) {
	//	sample = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
	//}
	//else if (rnd < w_diffuse + w_metal) {
	//	sample = operator()(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
 //   }
	//else if (rnd < w_diffuse + w_metal + w_clearcoat) {
	//	sample = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
	//}
	//else {
	//	sample = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
	//}

    if (rnd < w_glass) {
        sample = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    }
    else if (rnd < w_glass + w_diffuse) {
        sample = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    }
    else if (rnd < w_glass + w_diffuse + w_metal) {
        sample = operator()(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    }
    else {
        sample = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
    }

    return sample;
    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
