#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
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
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    // Check if incident inside (inside = only glass is needed, other weights all 0)
    bool isInside = dot(vertex.geometric_normal, dir_in) < 0;
    Spectrum f_glass = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    if (isInside) {
        return f_glass; 
    }

    // Compute remaining DisneyBSDF's components
    // Evaluate remaining key parameters for remaining components
	Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

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
	Spectrum F_m = C_0 + (1 - C_0) * pow(1 - h_dot_out, 5);
    // 2. D_m
    Real D_m = GTR2(h_local, roughness, anisotropic);
    // 3. G_m
    Real G_m = smith_masking_gtr2(to_local(frame, dir_in), roughness, anisotropic) *
               smith_masking_gtr2(to_local(frame, dir_out), roughness, anisotropic);
    // final f_metal
    f_metal = F_m * D_m * G_m / (4 * abs(n_dot_in));

    // Compute weight for each component
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_sheen = (1 - metallic) * sheen;
    Real w_metal = (1 - specular_transmission * (1 - metallic));
	Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;

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
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    // Check if incident inside (inside = only glass is needed, other weights all 0)
    bool isInside = dot(vertex.geometric_normal, dir_in) < 0;
    Real p_glass = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    if (isInside) {
        return p_glass;
    }

    // Compute remaining DisneyBSDF's components (DisneySheen is negelected in pdf)
    // Evaluate remaining key parameters for remaining components
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real p_diffuse = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    Real p_clearcoat = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
	Real p_metal = operator()(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic });

    // Compute weight for each component
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;

	// Compute the final result
    return w_diffuse * p_diffuse +
           w_metal * p_metal +
           w_clearcoat * p_clearcoat +
           w_glass * p_glass;
    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
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
    if (rnd < w_diffuse) {
		sample = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
	}
	else if (rnd < w_diffuse + w_metal) {
		sample = operator()(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    }
	else if (rnd < w_diffuse + w_metal + w_clearcoat) {
		sample = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
	}
	else {
		sample = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
	}

    return sample;
    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
