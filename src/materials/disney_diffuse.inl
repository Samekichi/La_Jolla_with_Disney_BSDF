Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
	
    // Evaluate key parameters for DisneyDiffuse
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    subsurface = std::clamp(subsurface, Real(0.01), Real(1));

    // Compute DisneyDiffuse's subcomponents
    // 0. common vars
    Vector3 h = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);  // cos(theta_in)
    Real n_dot_out = dot(frame.n, dir_out);  // cos(theta_out)
	Real h_dot_out = dot(h, dir_out);  // cos(theta_half_out)
    // 1. f_base_diffuse
    Real FD_90 = 0.5 + 2 * roughness * h_dot_out * h_dot_out;
    Real FD_in = 1 + (FD_90 - 1) * pow(1 - n_dot_in, 5);
    Real FD_out = 1 + (FD_90 - 1) * pow(1 - n_dot_out, 5);
    Spectrum f_base_diffuse = (base_color / c_PI) * (FD_in * FD_out * abs(n_dot_out));
    // 2. f_subsurface
    Real FSS_90 = roughness * pow(h_dot_out, 2);
    Real FSS_in = 1 + (FSS_90 - 1) * pow(1 - abs(n_dot_in), 5);
	Real FSS_out = 1 + (FSS_90 - 1) * pow(1 - abs(n_dot_out), 5);
    Spectrum f_subsurface = ((1.25 * base_color) / c_PI) * 
                            (FSS_in * FSS_out * (1 / (abs(n_dot_in) + abs(n_dot_out)) - 0.5) + 0.5) * 
                            n_dot_out;
    
	// Compute the final result
	Spectrum f_diffuse = (1 - subsurface) * f_base_diffuse + subsurface * f_subsurface;
    return f_diffuse;
    // return f_base_diffuse;
    // return f_subsurface;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    
    // Same as Lambertian, we importance sample the cosine hemisphere domain.
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
    // return Real(0);
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    
    // Same as Lambertian, we importance sample the cosine hemisphere domain.
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool) /* roughness */ };
    //return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
