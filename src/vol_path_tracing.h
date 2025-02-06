#pragma once

inline int update_medium(const Ray& ray, const PathVertex& vertex) {
    /*
        At medium transition. Update medium.
        - returns the new medium the `ray` will be in, right after passing through the surface at `vertex`
    */
    if (dot(ray.dir, vertex.geometric_normal) > 0) {
        return vertex.exterior_medium_id;
    }
    else {
        return vertex.interior_medium_id;
    }
}

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    /* From path_tracing.h: path_tracing()*/
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray camera_ray = sample_primary(scene.camera, screen_pos);
    // Ray differential for volumetric scattering is an unsolved problem,
    // so we disable it for volumetric path tracing for now.
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    std::optional<PathVertex> vertex_ = intersect(scene, camera_ray, ray_diff);
    // For HW 2's volumetric path tracing, we assume there is no environment map in the scene.
    if (!vertex_) {
        return make_zero_spectrum();
    }
    PathVertex vertex = *vertex_;
    /* end from path_tracing() */

    // Prepare related volumetric parameters at the intersection point
    Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], vertex.position);
    Real t = distance(camera_ray.org, vertex.position);
    // Compute the transmittance and emission at the intersection point
    Spectrum transmittance = exp(-sigma_a * t);
    Spectrum Le = make_zero_spectrum();
    if (is_light(scene.shapes[vertex.shape_id])) {
        Le = emission(vertex, -camera_ray.dir, scene);
    }
    return transmittance * Le;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    /* From path_tracing.h: path_tracing()*/
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray camera_ray = sample_primary(scene.camera, screen_pos);
    // Ray differential for volumetric scattering is an unsolved problem,
    // so we disable it for volumetric path tracing for now.
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    std::optional<PathVertex> vertex_ = intersect(scene, camera_ray, ray_diff);
    // For HW 2's volumetric path tracing, we assume there is no environment map in the scene.
    // no intersection -> still need to consider volumetric scattering -> cannot directly return 0 Spectrum
    PathVertex vertex = vertex_.value_or(PathVertex{});  
    /* end from path_tracing() */

    // Prepare related volumetric parameters at the intersection point
    Real u = next_pcg32_real<Real>(rng);
    Spectrum sigma_s = get_sigma_s(scene.media[scene.camera.medium_id], vertex.position);
    Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], vertex.position);
    Spectrum sigma_t_vec = sigma_s + sigma_a;
    Real sigma_t = sigma_t_vec.x;  // for this task, we assume sigma_t is monochromatic
    Real t = -log(1 - u) / sigma_t;
    Real t_hit = vertex_ ? distance(camera_ray.org, vertex.position) : infinity<Real>();  // if no intersection, t_hit is infinity
    
    if (t < t_hit) {
        // The sampled t is inside the volume,
        // account for volume scattering.
        Real trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t_vec * t);

#       // Compute L_scatter1(p, w) using Monte Carlo sampling,
        // i.e., Randomly sample 1 light source and 1 point on that light source's surface.
        Vector3 p = camera_ray.org + t * camera_ray.dir;
        
        // sample 1 light source
        Real light_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light& light = scene.lights[light_id];
        // sample 1 point on that light source's surface
        Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
        Real shape_w = next_pcg32_real<Real>(rng);
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);

        // compute L_scatter1(p, w)
        Real L_scatter1_pdf = light_pmf(scene, light_id) *
            pdf_point_on_light(light, point_on_light, p, scene);  // P(sample the light) * P(sample the point on the light)
        // i. compute phase function
        Vector3 dir_light = normalize(point_on_light.position - p);
        Spectrum phase = eval(get_phase_function(scene.media[scene.camera.medium_id]), camera_ray.dir, dir_light);
        // ii. compute Le at light
        Spectrum Le = emission(light, -dir_light, ray_diff.spread, point_on_light, scene);
        // *. visible check
        Ray shadow_ray{ p, dir_light,
                        get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                            distance(point_on_light.position, p) };
        if (occluded(scene, shadow_ray)) {
            return make_zero_spectrum();
        }
        // iii. compute L_scatter1(p, w)
        Spectrum L_scatter1_estimate = phase * Le * exp(-sigma_t * distance(point_on_light.position, p)) * abs(dot(dir_light, point_on_light.normal)) / distance_squared(point_on_light.position, p);  // ignoring visible test

        return transmittance / trans_pdf * sigma_s * L_scatter1_estimate / L_scatter1_pdf;
    }
    else {
        // The sampled t hits a surface,
        // account for surface emission.
        Real trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t_vec * t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -camera_ray.dir, scene);
        }
        return (transmittance / trans_pdf) * Le;
    }
}

inline Spectrum next_event_estimation(Ray ray,  // incoming ray from previous node
    int current_medium,
    int bounces,
    const Scene& scene,
    pcg32_state& rng) {
    Vector3 original_p = ray.org;
    Vector3 p = original_p;
    // Sample a point on a light source for NEE
    // sample 1 light source
    Real light_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light& light = scene.lights[light_id];
    // sample 1 point on that light source's surface
    Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
    Real shape_w = next_pcg32_real<Real>(rng);
    PointAndNormal p_prime =
        sample_point_on_light(light, original_p, light_uv, shape_w, scene);

    // Compute transmittance to light. Skip through index-matching shapes.
    Vector3 dir_light = normalize(p_prime.position - original_p);
    Real T_light = 1;
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir = 1;  // for multiple importance sampling
    while (true) {
        Ray shadow_ray = Ray{ p,
                              normalize(p_prime.position - p),
                              get_shadow_epsilon(scene),
                              (1 - get_shadow_epsilon(scene)) * distance(p_prime.position, p) };
        // Find the next intersection point
        RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
        std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, ray_diff);
        PathVertex vertex = vertex_.value_or(PathVertex{});
        Real next_t = distance(p, p_prime.position);
        if (vertex_) {
            next_t = distance(p, vertex.position);
        }
        // Account for the transmittance to next_t
        Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium], p);
        Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium], p);
        Spectrum sigma_t_vec = sigma_s + sigma_a;
        Real sigma_t = sigma_t_vec.x;  // for this task, we assume sigma_t is monochromatic
        if (shadow_medium >= 0) {
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!vertex_) {
            // Nothing is blocking, we're done (i.e., we hit the light source)
            break;
        }
        else {
            // Something is blocking: is it an opaque sureface OR index-matching surface?
            if (vertex.material_id >= 0) {
                // Opaque surface: we're blocked
                return make_zero_spectrum();
            }
            // otherwise, it's Index-matching surface
            // we want to pass throuh -- this introduces one extra connection vertex
            shadow_bounces += 1;

            if (scene.options.max_depth != -1 && bounces + shadow_bounces >= scene.options.max_depth) {
                // We've reached the max depth
                return make_zero_spectrum();
            }

            shadow_medium = update_medium(shadow_ray, vertex);
            p = p + (next_t + get_intersection_epsilon(scene)) * shadow_ray.dir;
        }
    }

    if (T_light > 0) {
        // Compute contrib = T_light * G * rho * L / pdf_NEE
        // G
        Real G = max(-dot(dir_light, p_prime.normal), Real(0)) / distance_squared(p_prime.position, original_p);
        // rho
        PhaseFunction phase = get_phase_function(scene.media[current_medium]);
        Spectrum rho = eval(phase, ray.dir, dir_light);
        // L
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
        // pdf_NEE
        Real pdf_NEE = light_pmf(scene, light_id) *
            pdf_point_on_light(light, p_prime, original_p, scene);
        // -> contrib
        Spectrum contrib = T_light * G * rho * L / pdf_NEE;
        // Multiple importance sampling :
        // - It's also possible that a phase function sampling + multiple exponential sampling
        //   will reach the light source.
        // - We also need to multiply with G to convert phase function PDF to area measure.
        // Compute w = pdf_NEE^2 / (pdf_NEE^2 + pdf_phase^2)
        Real pdf_phase = pdf_sample_phase(phase, ray.dir, dir_light);
        Real w = pdf_NEE * pdf_NEE / (pdf_NEE * pdf_NEE + pdf_phase * pdf_phase);
        return w * contrib;
    }
    return make_zero_spectrum();
}


// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    // Sample a camera ray
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    // Variables for recursively computing L_scatter(p, w)
    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1);  // contrib(path) / p(path)
    Spectrum radiance = make_zero_spectrum();
    // - path has at most (scene.options.max_depth + 1) nodes;
    // - if max_depth =-1, then we can bounce infinite amount of time;
    // - max_depth = 2 corresponds to the 1-scattering case in `vol_path_tracing_2` (i.e., camera -> p(t) -> light)
    int bounces = 0;  
    // Ray differential for volumetric scattering is an unsolved problem,
    // so we disable it for volumetric path tracing for now.
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    while (true) {
        bool scatter = false;
        // Find the next intersection point
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = vertex_.value_or(PathVertex{});  // no intersection -> still need to consider volumetric scattering -> cannot directly return 0 Spectrum
        
        // Ray might not intersect a surface, but we might be in a volume
        Spectrum transmittance = make_const_spectrum(1);
        Real trans_pdf = 1;

        // Sample a distance `t` towards the intersection, and check if it hits the intersection
        if (current_medium >= 0) {
            // sample one step distance `t` s.t. p(t) ~ exp(-sigma_t * t)
            Real u = next_pcg32_real<Real>(rng);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_t_vec = sigma_s + sigma_a;
            Real sigma_t = sigma_t_vec.x;  // for this task, we assume sigma_t is monochromatic
            Real t = -log(1 - u) / sigma_t;
            Real t_hit = vertex_ ? distance(ray.org, vertex.position) : infinity<Real>();  // if no intersection, t_hit is infinity
            // if t < t_hit, ray scatter in the volume
            if (t < t_hit) {
                scatter = true;
                // compute transmittance and trans_pdf
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_vec * t);
                // update next iter's ray's origin (direction will be updated later)
                ray.org = ray.org + t * ray.dir;
            }
            // else t >= t_hit, ray hits a surface
            else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_vec * t_hit);
                ray.org = ray.org + (t_hit + get_intersection_epsilon(scene)) * ray.dir;
            }
            
        }
        // Update accumulated transmittance
        current_path_throughput *= (transmittance / trans_pdf);

        // Hits a surface: add its contributed emission Le
        if (!scatter) {
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }
            radiance += current_path_throughput * Le;
        }
        // Check if need to cut off the path (i.e., there's a max_depth and we reached it)
        if (scene.options.max_depth != -1 && bounces >= scene.options.max_depth) {
            break;
        }
        // Hit index-matching surface (volume boundaries): update current_medium & continue to next iteration
        if (!scatter ) {
            if (vertex.material_id == -1) {
                current_medium = update_medium(ray, vertex);
                bounces += 1;
                continue;
            }
        }

        // Sample next direction
        if (scatter) { 
            // sample next iter's direction based on phase_function
            Vector3 next_dir = sample_phase_function(get_phase_function(scene.media[current_medium]), ray.dir, Vector2{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) }).value_or(Vector3(0, 0, 0));
            PhaseFunction phase = get_phase_function(scene.media[current_medium]);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);

            //Spectrum nee_contrib = next_event_estimation(ray, current_medium, bounces, scene, rng);
            //radiance += current_path_throughput * nee_contrib;

            current_path_throughput *= eval(phase, ray.dir, next_dir) / pdf_sample_phase(phase, ray.dir, next_dir) * sigma_s;
            // update ray.dir
            ray.dir = next_dir;
        }
        else {
            break;  // Hit a surface -- we don't need to deal with this yet
        }

        // Russian Roulette!!!
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), Real(0.95));  // Follow path_trace()'s pattern in path_tracing.h
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }

        // update iter counter
        bounces += 1;
    } // end while loop for L_scatter(p, w)
    
    return radiance;
}


// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    /* !!!! Copied from vol_path_tracing_3() !!!! */
    // Sample a camera ray
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    // Variables for recursively computing L_scatter(p, w)
    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1);  // contrib(path) / p(path)
    Spectrum radiance = make_zero_spectrum();
    // - path has at most (scene.options.max_depth + 1) nodes;
    // - if max_depth =-1, then we can bounce infinite amount of time;
    // - max_depth = 2 corresponds to the 1-scattering case in `vol_path_tracing_2` (i.e., camera -> p(t) -> light)
    int bounces = 0;
    // Ray differential for volumetric scattering is an unsolved problem,
    // so we disable it for volumetric path tracing for now.
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    while (true) {
        bool scatter = false;
        // Find the next intersection point
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        PathVertex vertex = vertex_.value_or(PathVertex{});  // no intersection -> still need to consider volumetric scattering -> cannot directly return 0 Spectrum

        // Ray might not intersect a surface, but we might be in a volume
        Spectrum transmittance = make_const_spectrum(1);
        Real trans_pdf = 1;

        // Sample a distance `t` towards the intersection, and check if it hits the intersection
        if (current_medium >= 0) {
            // sample one step distance `t` s.t. p(t) ~ exp(-sigma_t * t)
            Real u = next_pcg32_real<Real>(rng);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_t_vec = sigma_s + sigma_a;
            Real sigma_t = sigma_t_vec.x;  // for this task, we assume sigma_t is monochromatic
            Real t = -log(1 - u) / sigma_t;
            Real t_hit = vertex_ ? distance(ray.org, vertex.position) : infinity<Real>();  // if no intersection, t_hit is infinity
            // if t < t_hit, ray scatter in the volume
            if (t < t_hit) {
                scatter = true;
                // compute transmittance and trans_pdf
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t_vec * t);
                // update next iter's ray's origin (direction will be updated later)
                ray.org = ray.org + t * ray.dir;
            }
            // else t >= t_hit, ray hits a surface
            else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t_vec * t_hit);
                ray.org = ray.org + (t_hit + get_intersection_epsilon(scene)) * ray.dir;
            }

        }
        // Update accumulated transmittance
        current_path_throughput *= (transmittance / trans_pdf);

        // Hits a surface: add its contributed emission Le
        if (!scatter) {
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }
            radiance += current_path_throughput * Le;
        }
        // Check if need to cut off the path (i.e., there's a max_depth and we reached it)
        if (scene.options.max_depth != -1 && bounces >= scene.options.max_depth) {
            break;
        }
        // Hit index-matching surface (volume boundaries): update current_medium & continue to next iteration
        if (!scatter) {
            if (vertex.material_id == -1) {
                current_medium = update_medium(ray, vertex);
                bounces += 1;
                continue;
            }
        }

        // Sample next direction
        if (scatter) {
            // sample next iter's direction based on phase_function
            Vector3 next_dir = sample_phase_function(get_phase_function(scene.media[current_medium]), ray.dir, Vector2{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) }).value_or(Vector3(0, 0, 0));
            PhaseFunction phase = get_phase_function(scene.media[current_medium]);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            
            Spectrum nee_contrib = next_event_estimation(ray, current_medium, bounces, scene, rng);
            radiance += current_path_throughput * nee_contrib;
            break;
            //current_path_throughput *= eval(phase, ray.dir, next_dir) / pdf_sample_phase(phase, ray.dir, next_dir) * sigma_s;

            //// update ray.dir
            //ray.dir = next_dir;
        }
        else {
            break;  // Hit a surface -- we don't need to deal with this yet
        }

        // Russian Roulette!!!
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), Real(0.95));  // Follow path_trace()'s pattern in path_tracing.h
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }

        // update iter counter
        bounces += 1;
    } // end while loop for L_scatter(p, w)

    return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
