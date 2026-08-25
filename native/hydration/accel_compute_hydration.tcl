# AUTO-GENERATED accelerated variant of compute_hydration (see gen_accel_compute_hydration.py).
# Phase A COG+projection and Phase C binning are replaced with calls to
# native/hydro_project; every other line is extracted VERBATIM from the
# real vmdhole.tcl (never edited). Requires $::HP_ACCEL_BIN and
# $::HP_ACCEL_BATCH_DIR to be set before calling.
proc ::VMDHole::_accel_compute_hydration {} {
    # Water-density + free-energy profile along the pore (Klesse, Rao, Sansom &
    # Tucker 2019, J Mol Biol 431:3353; Help > References). rho(s) = <waters in
    # bin>/(pi*R(s)^2*dz) per frame; E = -kT ln(rho/rho_bulk), both averaged
    # separately across frames.
    variable state
    variable results
    variable result_frames
    variable hydration_data
    set wsel [string trim $state(water_sel)]
    if {$wsel eq ""} { set wsel "water and oxygen" }
    # Bulk density: MEASURED from this trajectory, not assumed. Falls back to the
    # literature value only when the system has too little bulk water to measure.
    # See measure_bulk_density for why a constant is the wrong default.
    set bulk ""
    set _bulk_src "measured"
    set _bulk_mol -1
    catch {set _bulk_mol [resolve_molid]}
    if {$_bulk_mol >= 0 && [llength $result_frames] > 0} {
        # First few frames are enough: bulk density is a system property, not a
        # per-frame one, and this costs two atomselects per sampled water.
        set _bf [lrange $result_frames 0 4]
        if {[catch {measure_bulk_density $_bulk_mol $_bf $wsel} _bres]} {
            vmdcon -warn "VMDHole: bulk-density measurement failed ($_bres) - using the literature value."
            set bulk ""
        } else {
            set bulk $_bres
        }
    }
    if {![string is double -strict $bulk] || $bulk <= 0} {
        set bulk 0.0334; set _bulk_src "literature default - too little bulk water to measure"
    }
    set state(water_bulk) $bulk
    # Report the CANONICAL selection, not the user's: when no single-site variant
    # covers every residue the measurement silently falls back to the raw selection,
    # which is the one case where the answer still depends on how the selection was
    # written. Printing it is the only way that path is visible.
    set _bulk_csel $wsel
    if {$_bulk_mol >= 0} { catch {set _bulk_csel [_canonical_water_sel $_bulk_mol [lindex $result_frames 0] $wsel]} }
    vmdcon -info "VMDHole: bulk water density = [format %.5f $bulk] A^-3 ($_bulk_src; counted over \"$_bulk_csel\")"
    set dz [expr {[info exists state(water_dz)] ? $state(water_dz) : 1.0}]
    if {![string is double -strict $dz] || $dz <= 0} { set dz 1.0 }
    # Density sampling cap (Å): a local cylindrical probe radius. HOLE's pore radius
    # balloons to the endrad cap (~15 Å) at the bulk vestibules, so count/(π·R²·dz)
    # there assumes a huge pure-bulk cylinder and reads FALSELY 'dry' (the cylinder
    # pokes through the funnel walls). Capping the radius used for BOTH the water-
    # inclusion test and the bin volume turns density into a local axial probe:
    #   - true bulk  → ρ ≈ ρ_bulk  (probe sits in water)
    #   - wide vestibule → ρ ≈ ρ_bulk  (no longer dragged down by protein volume)
    #   - constriction (R < cap) → cap doesn't bind → true pore density (dewetting kept)
    # 0 disables the cap (legacy full-pore-radius behaviour).
    set dcap [expr {[info exists state(water_dens_cap)] ? $state(water_dens_cap) : 0.0}]
    # Default 0 = OFF. This caps the RADIUS used for the volume (despite the name),
    # so wherever the pore is wider than the cap, density is computed as if R = cap.
    # The old default of 5.0 A was active by default and had NO citation and no
    # derivation - the same failure as the removed 8x-bulk ceiling - while silently
    # altering density in every wide bin, which is most of the mouth region. It is
    # kept as an opt-in tool (a real concern: as the pore opens into a vestibule,
    # pi*R^2*dz grows fast and dilutes the cross-section average over space that is
    # no longer pore), but it is no longer applied unasked. CHAP does not correct for
    # large-R dilution at all, and CHAP mode already forced this to 0 - so defaulting
    # it off also makes the two modes agree.
    if {![string is double -strict $dcap] || $dcap < 0} { set dcap 0.0 }
    # Density probe FLOOR (Å): the symmetric MINIMUM counterpart to dcap above. At a
    # sub-Angstrom hydrophobic gate (real and expected - literature/this project's own
    # case study report ~0.36 Å constrictions), the bin
    # volume pi*R^2*dz becomes tiny, so rho = count/(N*vol) is EXTREMELY sensitive to a
    # single stray/rare water there - one incidental count can spike rho/rho_bulk past 1
    # even though the site is mostly dry, which is what makes such spikes "hard to
    # interpret" against the surrounding hydrophobic lining. Flooring R damps that
    # amplification without touching wide/bulk regions (where R is already >> floor).
    # 0 = off (legacy, matches all prior outputs exactly).
    set rfloor [expr {[info exists state(water_radius_floor)] ? $state(water_radius_floor) : 0.0}]
    if {![string is double -strict $rfloor] || $rfloor < 0} { set rfloor 0.0 }
    # water_poisson_floor (default ON = solid plugin default): applies the RULE-OF-THREE
    # density floor (3/(N·V), capped at bulk) to the per-frame energy, so a fully-dewetted
    # gate reports a CONSERVATIVE, sampling-limited barrier instead of a −ln(0) blow-up.
    # OFF = CHAP-faithful: CHAP's density is a strictly-positive KDE spline, so its −ln(ρ)
    # essentially never hits log(0) (its mendInfinities safety net almost never fires); this
    # emulate that with a small BOUNDED epsilon floor + no upper cap, NOT CHAP's literal
    # ±max_float (which would explode the frame-averaged energy). chap_mode sets this OFF.
    set pf_on [expr {![info exists state(water_poisson_floor)] || $state(water_poisson_floor)}]
    # water_energy_shift (default ON, editable in BOTH modes; CHAP also does this - it
    # shifts the profile so the mean of the two pore-mouth energies is zero, see
    # chap_trajectory_analysis.cpp:1765). OFF = raw −kT·ln(ρ) with CHAP's unit-dependent
    # baseline. Kept a standalone user knob (not part of the chap_mode bundle).
    set shift_on [expr {![info exists state(water_energy_shift)] || $state(water_energy_shift)}]
    # Gaussian KDE for the water density along the axis, ON by default (see this proc's
    # header for the bandwidth default and why). state(water_kde_bw) blank/"auto"/non-
    # numeric -> AMISE-optimal per frame; a positive number -> that fixed bandwidth.
    # Accumulating K_h(z-z_water)*dz into the bins makes the existing density formula
    # apply unchanged (the kernel integrates to 1/water).
    set use_kde [expr {[info exists state(water_kde)] && $state(water_kde)}]
    set bw_spec [expr {[info exists state(water_kde_bw)] ? [string trim $state(water_kde_bw)] : ""}]
    set bw_auto [expr {$bw_spec eq "" || [string equal -nocase $bw_spec "auto"] || \
        ![string is double -strict $bw_spec] || $bw_spec <= 0}]
    set bw 1.0
    if {!$bw_auto} { set bw $bw_spec }
    # Performance/robustness cap on the AMISE estimator's O(N^2) direct
    # density-derivative sum: a DETERMINISTIC, evenly-spaced subsample of at
    # most this many waters (sorted, fixed indices - see the thinning below;
    # no RNG is involved) is used ONLY to SOLVE for the bandwidth, and the
    # bandwidth alone. The actual KDE evaluation
    # below still uses every real water. Plugin-specific safety valve (like
    # water_dens_cap/water_radius_floor), not part of CHAP's own algorithm,
    # needed because CHAP's own "fast" O(N) approximate density-derivative
    # method (deliberately not ported - see _amise_phi_direct) is
    # what keeps THEIR per-frame cost bounded for large samples; this plugin's
    # direct method needs its own bound instead.
    set amise_cap 100
    # Channel axis for binning (see water_use_cpoint_axis in this proc's header). Off:
    # re-fits PCA independently per frame - a bin index then doesn't mean the same
    # physical location across frames on a fluctuating trajectory. On: the run's fixed
    # CPOINT/CVECT via _manifest_axis (matches HOLE's own "coord" column).
    set _want_cpaxis [expr {[info exists state(water_use_cpoint_axis)] && $state(water_use_cpoint_axis)}]
    set fixed_axis {}
    # axis_use selects HOW each frame's binning axis is derived:
    #   manifest_full - the run's ONE static CPOINT+CVECT for every frame (axis never drifts).
    #   manifest_dir  - the run's known CVECT DIRECTION for every frame (Track CPOINT drifts
    #                   only the ORIGIN, not the direction - see _manifest_cvect_drifts_per_
    #                   frame), with each frame's own centerline centroid as the origin. This
    #                   bins on HOLE's ACTUAL search direction instead of re-deriving a noisier
    #                   per-frame PCA estimate of a direction already known exactly.
    #   pca           - per-frame PCA of that frame's own centerline (no usable manifest, the
    #                   DIRECTION genuinely drifts (Stabilize CVECT/rotation), or the user
    #                   turned the CPOINT/CVECT anchor off).
    set axis_use "pca"
    # Carries a note about the anchor being stood down through to the FINAL
    # status message below (an intermediate set state(status) here would just
    # get overwritten once the actual computation finishes and reports its
    # own summary - vmdcon logs it immediately either way).
    set _anchor_note ""
    if {$_want_cpaxis} {
        set _odrift 0; set _cdrift 0
        foreach frame $result_frames {
            if {![dict exists $results $frame run_dir]} continue
            set _run_dir [dict get $results $frame run_dir]
            set fixed_axis [_manifest_axis $_run_dir]
            if {$fixed_axis ne {}} {
                set _odrift [_manifest_axis_drifts_per_frame $_run_dir]
                set _cdrift [_manifest_cvect_drifts_per_frame $_run_dir]
                break
            }
        }
        if {$fixed_axis eq {}} {
            set axis_use "pca"; set fixed_axis {}
            set _anchor_note " CPOINT/CVECT anchor requested but unavailable (no manifest) - used per-frame PCA instead."
            vmdcon -info "VMDHole: hydration -$_anchor_note"
        } elseif {!$_odrift && !$_cdrift} {
            # Static axis for the whole run: the one manifest CPOINT+CVECT lands every frame
            # on HOLE's exact "coord" grid.
            set axis_use "manifest_full"
        } elseif {!$_cdrift} {
            # Track CPOINT drifts the ORIGIN but NOT the direction, so the run's cvect is
            # still the true axis direction for every frame. Bin on it (HOLE's own search
            # direction) with each frame's centerline centroid as origin - strictly more
            # consistent than a per-frame PCA direction that wobbles with centerline noise.
            set axis_use "manifest_dir"
        } else {
            # The axis DIRECTION genuinely evolves per frame (Stabilize CVECT / rotation) -
            # the one static cvect no longer describes most frames, so per-frame PCA (fit
            # fresh from each frame's own spheres) is the correct choice here.
            set axis_use "pca"; set fixed_axis {}
            set _anchor_note " This run's axis DIRECTION moves per frame (Stabilize CVECT) - used per-frame PCA instead."
            vmdcon -info "VMDHole: hydration -$_anchor_note"
        }
    }
    set bincount [dict create]
    set binrsum  [dict create]
    set binrn    [dict create]
    set nfdata 0
    set nfwater 0
    set total_w 0
    set nframes [llength $result_frames]
    # Per-frame density storage: list of dicts, one per processed frame.
    # Each dict: {frame <f> bins <dict bin->count>}
    set perframe_raw {}
    # CHAP's anchor points (chap_trajectory_analysis.cpp:1592-1594): the most
    # extreme low/high axial position reached by ANY single frame's own
    # envelope (arcLengthLoSummary.min()/arcLengthHiSummary.max() in their
    # code) - fixed, global values used later to shift the whole energy
    # profile so it reads zero there.
    set anchor_lo_extreme 1e30
    set anchor_hi_extreme -1e30
    vmdcon -info "VMDHole: hydration starting — $nframes frame(s)  water=\"$wsel\"  bulk=$bulk /A^3  T=$state(water_temp) K"
    # PHASE A (sequential - the only part that needs VMD/atomselect): per
    # frame, find the qualifying waters' axial coordinates (qco). This is the
    # CHEAP ~15% of the total cost (measured); the expensive part is the
    # per-frame AMISE bandwidth solve, deferred to Phase B below so it can run
    # in parallel across frames - they are fully independent of each other.
    if {$bw_auto} {
        error "hydro_project accelerated path requires a FIXED water_kde_bw (not auto/AMISE) - see NOTES-hydration-accel.md"
    }
    set _hp_batch_dir $::HP_ACCEL_BATCH_DIR
    file mkdir $_hp_batch_dir
    set _hp_batch_lines {}
    set _hp_frame_order {}
    set fi 0
    set frame_qco_need_bw [dict create]
    set frame_ctx [dict create]
    set frame_range [dict create]   ;# per-frame centerline axial range [cmin,cmax] for the energy coverage gate
    set frame_radii [dict create]   ;# per-frame mean pore radius PER BIN - MUST be keyed by frame (see Phase C)
    set _wsel_validated 0
    foreach frame $result_frames {
        incr fi
        if {![dict exists $results $frame]} { continue }
        set rec [dict get $results $frame]
        if {![dict exists $rec sph_file]} { continue }
        set sph [dict get $rec sph_file]
        set centers [parse_sph_centerline $sph]
        if {[llength $centers] < 2} { continue }
        set molid [expr {[dict exists $rec draw_molid] ? [dict get $rec draw_molid] : [resolve_molid]}]
        # Validate the water selection ONCE, up front, on the first usable frame's
        # molid. A selection SYNTAX error is identical every frame; letting the
        # per-frame catch below swallow it and `continue` would count every frame
        # as usable-but-empty and report a completely DRY trajectory (false
        # dewetting) as a successful result. Fail loudly instead.
        if {!$_wsel_validated} {
            if {[catch {set _wt [atomselect $molid $wsel]} _werr]} {
                set hydration_data {}
                set state(status) "Hydration: water selection \"$wsel\" is invalid: $_werr"
                return 0
            }
            # RECOGNITION: a syntactically-valid selection that matches ZERO atoms means
            # there is no explicit water in this system (a bare PDB, or a dry / implicit-
            # solvent run) - distinct from a real dry pore. Say so and point at the
            # geometry-based estimates, instead of reporting a false all-dry result.
            set _wn 0; catch {set _wn [$_wt num]}
            catch {$_wt delete}
            if {$_wn == 0} {
                set hydration_data {}
                set state(status) "Hydration: no water matched \"$wsel\" in this system — it has no explicit water (a bare PDB or a dry/implicit run). The water-density / free-energy profile needs a solvated trajectory; use the geometric conductance/passability estimates (Passability on the Pore Profile tab) instead."
                return 0
            }
            set _wsel_validated 1
        }
        if {$axis_use eq "manifest_full"} {
            ;# _manifest_axis returns {cpx cpy cpz ux uy uz} - origin first, then axis.
            lassign $fixed_axis mx my mz ux uy uz
        } elseif {$axis_use eq "manifest_dir"} {
            ;# Direction = the run's known cvect (HOLE's search axis); origin = THIS frame's
            ;# centerline centroid (from oriented_axis) so Track's origin drift is absorbed.
            lassign [oriented_axis $centers] _pux _puy _puz mx my mz
            lassign $fixed_axis _ox _oy _oz ux uy uz
        } else {
            lassign [oriented_axis $centers] ux uy uz mx my mz
        }
        # envelope + channel bbox
        set env {}
        set mnx 1e20; set mxx -1e20; set mny 1e20; set mxy -1e20; set mnz 1e20; set mxz -1e20
        set maxr 0.0
        foreach c $centers {
            lassign $c cx cy cz r
            set co [expr {($cx-$mx)*$ux + ($cy-$my)*$uy + ($cz-$mz)*$uz}]
            lappend env [list $co $r]
            if {$cx < $mnx} {set mnx $cx}; if {$cx > $mxx} {set mxx $cx}
            if {$cy < $mny} {set mny $cy}; if {$cy > $mxy} {set mxy $cy}
            if {$cz < $mnz} {set mnz $cz}; if {$cz > $mxz} {set mxz $cz}
            if {$r > $maxr} {set maxr $r}
        }
        set env [lsort -real -index 0 $env]
        set cmin [lindex [lindex $env 0] 0]; set cmax [lindex [lindex $env end] 0]
        if {$cmin < $anchor_lo_extreme} { set anchor_lo_extreme $cmin }
        if {$cmax > $anchor_hi_extreme} { set anchor_hi_extreme $cmax }
        # per-frame envelope radius accumulation. binrsum/binrn build the
        # trajectory-MEAN radius per bin (used for the Density view, a valid linear
        # quantity). _frame_rsum/_frame_rn additionally build THIS frame's OWN mean
        # radius per bin, stored in perframe_raw below and used for the per-frame
        # ENERGY volume (C2) - a breathing/gated pore's energy needs each frame's
        # own R_t(z), not the trajectory mean (Jensen: mean(-ln(n/V_t)) != -ln(n/V̄)).
        set _frame_rsum [dict create]; set _frame_rn [dict create]
        foreach e $env {
            lassign $e co r
            set bi [expr {int(floor($co/$dz))}]
            dict set binrsum $bi [expr {([dict exists $binrsum $bi] ? [dict get $binrsum $bi] : 0.0) + $r}]
            dict incr binrn $bi
            dict set _frame_rsum $bi [expr {([dict exists $_frame_rsum $bi] ? [dict get $_frame_rsum $bi] : 0.0) + $r}]
            dict incr _frame_rn $bi
        }
        set _frame_radii [dict create]
        foreach _rb [dict keys $_frame_rsum] {
            dict set _frame_radii $_rb [expr {[dict get $_frame_rsum $_rb] / double([dict get $_frame_rn $_rb])}]
        }
        incr nfdata
        # waters: pre-filter to the channel bbox so the scan stays cheap
        set m [expr {$maxr + 2.0}]
        set q "($wsel) and x > [expr {$mnx-$m}] and x < [expr {$mxx+$m}] and y > [expr {$mny-$m}] and y < [expr {$mxy+$m}] and z > [expr {$mnz-$m}] and z < [expr {$mxz+$m}]"
        # The bare selection passed the up-front validation, so a failure HERE is a
        # real per-frame problem (e.g. the molecule was deleted mid-run) - surface
        # it instead of silently skipping the frame (which would under-count water).
        if {[catch {set wa [atomselect $molid $q frame $frame]} _werr]} {
            set hydration_data {}
            set state(status) "Hydration: water selection failed at frame $frame: $_werr"
            return 0
        }
        set wpos [$wa get {x y z}]
        set wres [$wa get residue]
        $wa delete
        if {[llength $wpos] > 0} { incr nfwater }
        set _hp_in [file join $_hp_batch_dir [format "f%06d.in" $frame]]
        set _hp_out [file join $_hp_batch_dir [format "f%06d.out" $frame]]
        set _hp_fh [open $_hp_in w]
        puts $_hp_fh [format "AXIS %.17g %.17g %.17g %.17g %.17g %.17g" $mx $my $mz $ux $uy $uz]
        puts $_hp_fh [format "RANGE %.17g %.17g" $cmin $cmax]
        puts $_hp_fh [format "DCAP %.17g" $dcap]
        set _hp_ecos {}; set _hp_ers {}
        foreach _hp_e $env { lassign $_hp_e _hp_co _hp_r; lappend _hp_ecos [format %.17g $_hp_co]; lappend _hp_ers [format %.17g $_hp_r] }
        puts $_hp_fh "ENV [llength $env]"
        puts $_hp_fh [join $_hp_ecos " "]
        puts $_hp_fh [join $_hp_ers " "]
        set _hp_xs {}; set _hp_ys {}; set _hp_zs {}
        foreach _hp_p $wpos { lassign $_hp_p _hp_x _hp_y _hp_z; lappend _hp_xs [format %.17g $_hp_x]; lappend _hp_ys [format %.17g $_hp_y]; lappend _hp_zs [format %.17g $_hp_z] }
        puts $_hp_fh "WATERS [llength $wpos]"
        puts $_hp_fh [join $wres " "]
        puts $_hp_fh [join $_hp_xs " "]
        puts $_hp_fh [join $_hp_ys " "]
        puts $_hp_fh [join $_hp_zs " "]
        puts $_hp_fh [format "BIN %.17g %.17g %d" $dz $bw $use_kde]
        close $_hp_fh
        lappend _hp_batch_lines "$_hp_in\t$_hp_out"
        lappend _hp_frame_order $frame
        dict set frame_ctx $frame {}
        dict set frame_range $frame [list $cmin $cmax]
        dict set frame_radii $frame $_frame_radii
        if {$fi % 20 == 0 || $fi == $nframes} {
            set state(status) "Hydration: scanning frame $fi / $nframes ..."
            update
        }
    }

    # PHASE B (parallel): solve the AMISE-optimal bandwidth for every frame
    # that needs one, across a Thread worker pool (falls back to serial
    # automatically for small frame counts or if Thread is unavailable).
    set frame_bw [dict create]
    if {[dict size $frame_qco_need_bw] > 0} {
        set state(status) "Hydration: solving bandwidths for [dict size $frame_qco_need_bw] frame(s) ([resolve_job_count] parallel)..."
        update
        set frame_bw [_amise_bandwidths_parallel $frame_qco_need_bw]
    }

    # PHASE C (sequential, cheap - pure arithmetic on already-collected qco):
    # bin every frame's waters using its own (now known) bandwidth.

    set frame_bw [dict create]
    set bincount [dict create]
    set _hp_batch_file [file join $_hp_batch_dir "batch.txt"]
    set _hp_global_file [file join $_hp_batch_dir "global.out"]
    set _hp_bf [open $_hp_batch_file w]
    foreach _l $_hp_batch_lines { puts $_hp_bf $_l }
    close $_hp_bf
    if {[llength $_hp_frame_order] > 0} {
        set _hp_cmd "[shell_quote $::HP_ACCEL_BIN] --bin --bin-global [shell_quote $_hp_global_file] --batch [shell_quote $_hp_batch_file]"
        if {[catch {exec sh -c $_hp_cmd} _hp_err]} {
            error "hydro_project batch failed: $_hp_err"
        }
        set _hp_gf [open $_hp_global_file r]
        set _hp_gline1 [gets $_hp_gf]
        set _hp_gline2 [gets $_hp_gf]
        close $_hp_gf
        if {[lindex $_hp_gline1 0] ne "BINS"} { error "hydro_project: malformed global output" }
        set _hp_glo [lindex $_hp_gline1 1]; set _hp_ghi [lindex $_hp_gline1 2]
        if {$_hp_ghi >= $_hp_glo} {
            set _hp_gvals $_hp_gline2
            set _hp_bi $_hp_glo
            foreach _hp_v $_hp_gvals { dict set bincount $_hp_bi $_hp_v; incr _hp_bi }
        }
        foreach frame $_hp_frame_order {
            set _hp_out [file join $_hp_batch_dir [format "f%06d.out" $frame]]
            set _hp_fh [open $_hp_out r]
            set _hp_l1 [gets $_hp_fh]
            set _hp_l2 [gets $_hp_fh]
            set _hp_l3 [gets $_hp_fh]
            set _hp_l4 {}
            if {[lindex $_hp_l3 0] eq "BINS"} { set _hp_l4 [gets $_hp_fh] }
            close $_hp_fh
            incr total_w [lindex $_hp_l1 1]
            set frame_bincount [dict create]
            if {[lindex $_hp_l3 0] eq "BINS"} {
                set _hp_lo [lindex $_hp_l3 1]; set _hp_hi [lindex $_hp_l3 2]
                if {$_hp_hi >= $_hp_lo} {
                    set _hp_bi $_hp_lo
                    foreach _hp_v $_hp_l4 { dict set frame_bincount $_hp_bi $_hp_v; incr _hp_bi }
                }
            }
            set _fcmin -1e30; set _fcmax 1e30
            if {[dict exists $frame_range $frame]} { lassign [dict get $frame_range $frame] _fcmin _fcmax }
            set _fr [expr {[dict exists $frame_radii $frame] ? [dict get $frame_radii $frame] : [dict create]}]
            lappend perframe_raw [dict create frame $frame bins $frame_bincount radii $_fr cmin $_fcmin cmax $_fcmax]
        }
    }
    set _T $state(water_temp)
    if {![string is double -strict $_T] || $_T <= 0} { set _T 300 }
    set kT [expr {1.9872041e-3 * $_T}]
    set coords {}; set occ {}; set dens {}; set countspf {}
    set meanrs {}
    set clamp_los {}
    set sorted_bins_all [lsort -integer [dict keys $binrn]]
    foreach bi $sorted_bins_all {
        set nr [dict get $binrn $bi]
        if {$nr == 0} { continue }
        set meanr [expr {[dict get $binrsum $bi] / double($nr)}]
        set meanr_eff [expr {($dcap > 0 && $meanr > $dcap) ? $dcap : $meanr}]
        if {$rfloor > 0 && $meanr_eff < $rfloor} { set meanr_eff $rfloor }
        set vol   [expr {3.14159265 * $meanr_eff * $meanr_eff * $dz}]
        set cnt   [expr {[dict exists $bincount $bi] ? [dict get $bincount $bi] : 0}]
        set perfr [expr {$cnt / double($nfdata)}]
        # NOTE: occupancy/density are built per-frame below (perframe_raw), not as
        # <N>/(pi <R>^2 dz) from the mean radius here. $vol/$meanr survive only for the
        # per-bin Poisson clamp and the reported mean-radius profile, which ARE per-bin
        # quantities.
        # Energy density floor = the RULE OF THREE (Hanley & Lippman-Hand, JAMA 1983;
        # 249:1743): observing ZERO waters in N frames of bin volume V bounds the true
        # density at ≈3/(N·V) at 95% confidence, capping the free-energy barrier at
        # −kT·ln(3/(N·V)) - a citable convention preventing −ln(0), reporting a
        # fully-dewetted gate as a conservative lower bound. In absolute Å⁻³ to match
        # CHAP's raw −ln(density) energy.
        # CAPPED at bulk: without this, a narrow/under-sampled bin's 3/(N·V) can exceed
        # bulk, flooring a dry gate above bulk and, after the mouth-zero shift, reading
        # spuriously WET. min(...,bulk) makes an under-sampled gate read ~0 instead.
        # pf_on (default): rule-of-three floor capped at bulk. pf OFF (chap_mode): a small
        # bounded epsilon so −ln stays finite, no bulk cap.
        set _clamp_lo [expr {$pf_on \
            ? (($nfdata > 0 && $vol > 0) ? min(3.0 / ($nfdata * $vol), $bulk) : 0.0065 * $bulk) \
            : 1.0e-4 * $bulk}]
        lappend coords  [expr {($bi + 0.5) * $dz}]
        lappend countspf $perfr
        lappend clamp_los $_clamp_lo
        lappend meanrs  $meanr
    }
    # $energy is built from a PER-FRAME G (below), not a log of the already-averaged
    # $occ above.
    #
    # PERFORMANCE (10-30k frame trajectories): a single pass over perframe_raw with a
    # running (Welford) mean/variance accumulator per bin, not a dense frames x bins
    # matrix (which would cost tens-hundreds of MB for no benefit - every value is only
    # needed once). Every frame must still touch every bin (skipping "didn't reach this
    # bin" would bias the mean upward), so time stays O(frames x bins) either way.
    set sorted_bins $sorted_bins_all
    set nb [llength $coords]
    set bi_lo_anchor [expr {int(floor($anchor_lo_extreme/$dz))}]
    set bi_hi_anchor [expr {int(floor($anchor_hi_extreme/$dz))}]
    set idx_lo_anchor -1; set idx_hi_anchor -1
    for {set bi2 0} {$bi2 < $nb} {incr bi2} {
        set bi_actual [lindex $sorted_bins $bi2]
        if {$bi_actual == $bi_lo_anchor} { set idx_lo_anchor $bi2 }
        if {$bi_actual == $bi_hi_anchor} { set idx_hi_anchor $bi2 }
    }
    # Precompute each bin's volume once (was recomputed per-frame before).
    set bin_vol {}
    foreach bi $sorted_bins {
        set nr [dict get $binrn $bi]
        if {$nr == 0} { lappend bin_vol 0.0; continue }
        set meanr [expr {[dict get $binrsum $bi] / double($nr)}]
        set meanr_eff [expr {($dcap > 0 && $meanr > $dcap) ? $dcap : $meanr}]
        if {$rfloor > 0 && $meanr_eff < $rfloor} { set meanr_eff $rfloor }
        lappend bin_vol [expr {3.14159265 * $meanr_eff * $meanr_eff * $dz}]
    }
    # Running Welford accumulators, O(bins) memory regardless of frame count.
    set occ_n [lrepeat $nb 0]; set occ_mean [lrepeat $nb 0.0]; set occ_m2 [lrepeat $nb 0.0]
    set g_n   [lrepeat $nb 0]; set g_mean   [lrepeat $nb 0.0]; set g_m2   [lrepeat $nb 0.0]
    set anchor_lo_gvals {}
    set anchor_hi_gvals {}
    # perframe_occ/perframe_frames (the dense frames x bins matrix) is built
    # whenever the Hydration tab is in the per-frame heatmap view OR the
    # trajectory is small enough that the memory cost is negligible either
    # way - a view-only gate would mean Compute-while-viewing-Density then
    # switching to Per-frame afterward shows a permanent "not available"
    # placeholder, even on a trajectory small enough that the O(frames x bins)
    # memory this gate exists to avoid is negligible. That cost only matters
    # at the 10-30k-frame scale, so the size threshold below is chosen
    # generously above ordinary trajectory lengths while staying well under it.
    set want_pf_matrix [expr {([info exists state(hydration_view)] && $state(hydration_view) in {heatmap heatmap_g}) \
        || $nframes <= 2000}]
    set perframe_occ {}
    set perframe_frames {}
    foreach pf_entry $perframe_raw {
        set pf_bins [dict get $pf_entry bins]
        set pf_radii [expr {[dict exists $pf_entry radii] ? [dict get $pf_entry radii] : [dict create]}]
        set _fcmin [expr {[dict exists $pf_entry cmin] ? [dict get $pf_entry cmin] : -1e30}]
        set _fcmax [expr {[dict exists $pf_entry cmax] ? [dict get $pf_entry cmax] : 1e30}]
        if {$want_pf_matrix} { set pf_occ {} }
        for {set bi2 0} {$bi2 < $nb} {incr bi2} {
            set bi   [lindex $sorted_bins $bi2]
            set vol  [lindex $bin_vol $bi2]
            set cnt  [expr {[dict exists $pf_bins $bi] ? [dict get $pf_bins $bi] : 0}]
            # PER-FRAME VOLUME for the density too (R-03). CHAP computes the density
            # once PER FRAME from that frame's own R(s) and then averages the per-frame
            # profiles - stated in the paper (Fig. 4 caption: "carried out individually
            # for each frame ... averaged over all time steps") and done in the source
            # (chap_trajectory_analysis.cpp: NumberDensityCalculator inside the per-frame
            # loop, then SummaryStatistics::mean()). Dividing the MEAN count by the MEAN
            # radius' volume is a different quantity, because 1/R^2 is convex and the pore
            # breathes (Jensen again - the same trap as the per-frame energy). Measured on
            # a 50-frame GABA run: ~2% median, 8-9% at the constriction, up to 34% in the
            # most strongly breathing bin. The fallback below is unreachable for any bin
            # inside [cmin,cmax] (HOLE samples at 0.5 A against dz=1.0, so every in-range
            # bin holds at least one centerline point - measured 0 gaps in 6338 bin-frames);
            # it only covers out-of-range bins, which the coverage gate then drops anyway.
            if {[dict exists $pf_radii $bi]} {
                set _fr [dict get $pf_radii $bi]
                set _fr [expr {($dcap > 0 && $_fr > $dcap) ? $dcap : $_fr}]
                if {$rfloor > 0 && $_fr < $rfloor} { set _fr $rfloor }
                set _vol_e [expr {3.14159265 * $_fr * $_fr * $dz}]
            } else {
                set _vol_e $vol
            }
            set rho  [expr {($_vol_e > 0) ? $cnt / $_vol_e : 0.0}]
            set v    [expr {($bulk > 0) ? $rho / $bulk : 0.0}]
            if {$want_pf_matrix} { lappend pf_occ $v }
            # ENERGY/DENSITY COVERAGE GATE: a frame contributes to a bin's statistics only if
            # the bin lies within THAT frame's own centerline range [cmin,cmax] - a frame whose
            # centerline didn't reach this bin has NO measurement here, not a dry zero. Matches
            # CHAP, which averages each bin only over frames whose own pathway reaches it. The
            # PER-FRAME HEATMAP above is deliberately ungated (shows each frame's own value).
            set _co [expr {($bi + 0.5) * $dz}]
            if {$_co < $_fcmin || $_co > $_fcmax} { continue }
            # Energy shares the density's per-frame volume computed above: both are
            # -kT ln / a division by pi R_t^2 dz with R_t the FRAME's own radius, so they
            # must use the identical convention or THIS FRAME's gv would not be the
            # Boltzmann inverse of THIS FRAME's own density. (The two are still averaged
            # separately across frames afterward - see the note below.)
            set rho_e $rho
            # Energy from RAW per-frame density (CHAP's own -ln(density); no bulk
            # division). Clamps in the SAME absolute (Å⁻³) units: $clo (from
            # clamp_los) is a raw-density Poisson floor, upper bound 8·bulk.
            set rc $rho_e
            set clo [lindex $clamp_los $bi2]
            # NO upper cap, in either mode - a wide-open pore's raw density is real geometry,
            # not an artifact to clamp (CHAP doesn't clamp it either; its answer is the KDE
            # bandwidth plus the red-hatched R<1.4A annotation this plugin already draws).
            # The 1e6*bulk below is an overflow guard for log(), not a physical bound.
            set chi [expr {1.0e6 * $bulk}]
            if {$rc > $chi} { set rc $chi }
            # Both modes floor at clo before taking the log - $clo is already the
            # right value for whichever mode this is (rule-of-three for pf_on, a
            # small bounded epsilon for CHAP mode; see clo's own definition above).
            # Must NOT reproduce CHAP's literal mendInfinities() fallback
            # (+/- numeric_limits<real>::max() for an exactly-zero-density
            # frame): that relies on CHAP's KDE-smoothed density essentially
            # never hitting exact zero, but this plugin's per-frame bin COUNT
            # is not a continuous KDE spline - a genuinely dry frame at a
            # genuinely dry gate (the case this tool exists to find) can hit
            # an exact zero count on real trajectories, and a single such
            # frame contributing +3.4e38 kcal/mol into a Welford mean would
            # explode that bin's averaged energy and dwarf every other bin on
            # the plot. The bounded-epsilon floor (see clo's own comment)
            # avoids that.
            if {$rc < $clo} { set rc $clo }
            set gv [expr {-1.0 * $kT * log($rc)}]

            set n [expr {[lindex $occ_n $bi2] + 1}]
            set delta [expr {$v - [lindex $occ_mean $bi2]}]
            set mean [expr {[lindex $occ_mean $bi2] + $delta/double($n)}]
            set m2 [expr {[lindex $occ_m2 $bi2] + $delta*($v-$mean)}]
            lset occ_n $bi2 $n; lset occ_mean $bi2 $mean; lset occ_m2 $bi2 $m2

            set gn [expr {[lindex $g_n $bi2] + 1}]
            set gdelta [expr {$gv - [lindex $g_mean $bi2]}]
            set gmean [expr {[lindex $g_mean $bi2] + $gdelta/double($gn)}]
            set gm2 [expr {[lindex $g_m2 $bi2] + $gdelta*($gv-$gmean)}]
            lset g_n $bi2 $gn; lset g_mean $bi2 $gmean; lset g_m2 $bi2 $gm2

            if {$bi2 == $idx_lo_anchor} { lappend anchor_lo_gvals $gv }
            if {$bi2 == $idx_hi_anchor} { lappend anchor_hi_gvals $gv }
        }
        if {$want_pf_matrix} {
            lappend perframe_occ $pf_occ
            lappend perframe_frames [dict get $pf_entry frame]
        }
    }
    # Energy is mean(-kT ln rho_frame) (g_mean below), not -kT ln(mean rho) - the same
    # per-frame convention CHAP uses (BoltzmannEnergyCalculator inside its own
    # per-frame loop; chap_trajectory_analysis.cpp).
    set energy {}; set occ_std {}; set energy_std {}
    set occ {}; set dens {}
    for {set bi2 0} {$bi2 < $nb} {incr bi2} {
        set _om [lindex $occ_mean $bi2]
        lappend occ  $_om
        lappend dens [expr {$_om * $bulk}]
        lappend energy [lindex $g_mean $bi2]
        set nv [lindex $occ_n $bi2]
        lappend occ_std    [expr {$nv > 1 ? sqrt([lindex $occ_m2 $bi2] / double($nv-1)) : 0.0}]
        lappend energy_std [expr {$nv > 1 ? sqrt([lindex $g_m2   $bi2] / double($nv-1)) : 0.0}]
    }
    # VACUOUS-BIN TRIM: a bin where NO frame's own HOLE centerline ever reached it
    # (occ_n==0) has occ_mean/g_mean stuck at their Welford INITIAL value 0.0 - not
    # a measurement. Reporting that as occ=0 ("totally dry") AND G=0 ("exactly
    # bulk") for the same non-measurement is self-contradictory (the cliff-to-zero
    # seen at both ends of a real energy plot). Distinct from a real bin sitting at
    # the rule-of-three floor (std=0, occ_n>0, a genuine citable measurement) -
    # this trim removes only zero-coverage bins, which by construction occur only
    # at the geometric extremes (see the R-03 coverage-gate note above). Not the
    # same class of change as the y-range trim below: this removes a
    # non-measurement, that never removes real data.
    set _keep {}
    for {set bi2 0} {$bi2 < $nb} {incr bi2} {
        if {[lindex $occ_n $bi2] > 0} { lappend _keep $bi2 }
    }
    if {[llength $_keep] < $nb} {
        set _nc {}; set _no {}; set _nd {}; set _ne {}; set _nos {}; set _nes {}
        set _ncp {}; set _nmr {}
        foreach bi2 $_keep {
            lappend _nc  [lindex $coords $bi2]
            lappend _no  [lindex $occ $bi2]
            lappend _nd  [lindex $dens $bi2]
            lappend _ne  [lindex $energy $bi2]
            lappend _nos [lindex $occ_std $bi2]
            lappend _nes [lindex $energy_std $bi2]
            lappend _ncp [lindex $countspf $bi2]
            lappend _nmr [lindex $meanrs $bi2]
        }
        set coords $_nc; set occ $_no; set dens $_nd; set energy $_ne
        set occ_std $_nos; set energy_std $_nes; set countspf $_ncp; set meanrs $_nmr
        # perframe_occ rows are built over the SAME bi2 range (see the per-frame
        # loop above) - keep them index-aligned with the now-trimmed coords, or
        # every pfdens/Per-frame-heatmap consumer that pairs coords[i] with
        # perframe_occ[frame][i] by index would silently misalign.
        if {$want_pf_matrix} {
            set _npo {}
            foreach _row $perframe_occ {
                set _nr {}
                foreach bi2 $_keep { lappend _nr [lindex $_row $bi2] }
                lappend _npo $_nr
            }
            set perframe_occ $_npo
        }
        set nb [llength $coords]
    }
    # Apply the shift (a constant added to every bin's mean energy does not
    # change energy_std - Var(X+c) = Var(X)).
    set shift 0.0
    if {$shift_on && [llength $anchor_lo_gvals] > 0 && [llength $anchor_hi_gvals] > 0} {
        set alo 0.0
        foreach v $anchor_lo_gvals { set alo [expr {$alo + $v}] }
        set alo [expr {$alo / double([llength $anchor_lo_gvals])}]
        set ahi 0.0
        foreach v $anchor_hi_gvals { set ahi [expr {$ahi + $v}] }
        set ahi [expr {$ahi / double([llength $anchor_hi_gvals])}]
        set shift [expr {-0.5 * ($alo + $ahi)}]
        set _shifted {}
        foreach g $energy { lappend _shifted [expr {$g + $shift}] }
        set energy $_shifted
    }
    set _gmin 0.0; set _gmax 0.0
    foreach _g $energy {
        if {$_g < $_gmin} { set _gmin $_g }
        if {$_g > $_gmax} { set _gmax $_g }
    }
    # axis_mode/axis record which channel axis these coords are registered to, so
    # gz_profile/dens_profile (which project a DISPLAY frame's own centerline onto
    # this SAME coordinate system to overlay colours) stay consistent with how the
    # bins were built rather than silently re-fitting their own per-frame PCA axis.
    #   cpoint = static manifest CPOINT+CVECT; cvect = manifest DIRECTION + per-frame origin
    #   (NOTE: for "cvect" the stored `axis` carries the manifest ORIGIN, which binning did
    #    NOT use - _hydration_display_axis deliberately reads only the DIRECTION from it and
    #    recomputes the per-frame centroid origin. Do not "fix" it to use the stored origin,
    #    that would break registration); pca = per-frame PCA.
    set _axis_mode [dict get {manifest_full cpoint manifest_dir cvect pca pca} $axis_use]
    set hydration_data [dict create coords $coords occupancy $occ energy $energy \
        occ_std $occ_std energy_std $energy_std \
        density $dens countspf $countspf radii $meanrs nframes $nfdata nframes_water $nfwater \
        total_waters $total_w wsel $wsel bulk $bulk dz $dz kT $kT \
        min_g $_gmin max_g $_gmax energy_shift $shift \
        axis_mode $_axis_mode axis $fixed_axis \
        perframe_occ $perframe_occ perframe_frames $perframe_frames]
    # Persist the computed profile next to the results so it can be reloaded on a later
    # import without recomputing (load_hydration_for_root). The dict is a plain string.
    catch {
        set _f0 [lindex $result_frames 0]
        set _root [file dirname [dict get $results $_f0 run_dir]]
        set _hf [open [file join $_root vmdhole_hydration.dat] w]
        puts $_hf $hydration_data
        close $_hf
    }
    # Invalidate gz/dens/pfdens-scheme colored surfaces so they rebuild with new
    # hydration_data. pfdens depends on perframe_occ specifically - a stale cached
    # surface here would keep showing a PREVIOUS run's per-frame water pattern
    # after a recompute changed the underlying occupancy, with no sign anything
    # was wrong (the file's mtime alone can't tell hydration_data changed).
    variable plot_cache
    variable plot_cache_order
    if {[info exists results]} {
        foreach _f [dict keys $results] {
            set _rd [dict get [dict get $results $_f] run_dir]
            catch {foreach _p [glob -nocomplain \
                    [file join $_rd "hole_hydro_gz*.vmd_plot"] \
                    [file join $_rd "hole_hydro_dens*.vmd_plot"] \
                    [file join $_rd "hole_hydro_pfdens*.vmd_plot"] \
                    [file join $_rd "hole_hydro_values_gz*.dat"] \
                    [file join $_rd "hole_hydro_values_dens*.dat"] \
                    [file join $_rd "hole_hydro_values_pfdens*.dat"]] {
                file delete $_p
            }}
        }
    }
    # Also invalidate the MEAN PROFILE's gz/dens-coloured surface files: their filename tag
    # keys on plot_data_version (bumped only by a new HOLE RUN), NOT on hydration, so a
    # hydration RECOMPUTE with no re-run would otherwise serve the OLD hydration colouring
    # from disk, so the mean profile can't show old hydration. Delete both the
    # smoothed display files (_prop_gz/dens) and the smooth-independent cache (_propus_gz/
    # dens); the mean surface then rebuilds fresh below. mean_dir is the same for every frame
    # (parent of run_dir), so derive it once.
    if {[info exists results] && [dict size $results] > 0} {
        set _mf0 [lindex [dict keys $results] 0]
        catch {
            set _mdir [file join [file dirname [dict get [dict get $results $_mf0] run_dir]] mean_profile]
            foreach _p [glob -nocomplain \
                    [file join $_mdir "mean_profile_*_prop_gz*.vmd_plot"] \
                    [file join $_mdir "mean_profile_*_prop_dens*.vmd_plot"] \
                    [file join $_mdir "mean_profile_*_propus_gz*.vmd_plot"] \
                    [file join $_mdir "mean_profile_*_propus_dens*.vmd_plot"]] {
                file delete $_p
            }
        }
    }
    foreach _k [array names plot_cache] {
        if {[string match "*hole_hydro_gz*" $_k] || [string match "*hole_hydro_dens*" $_k] || \
                [string match "*hole_hydro_pfdens*" $_k]} {
            catch {unset plot_cache($_k)}
        }
    }
    set plot_cache_order [lsearch -all -not -inline -regexp $plot_cache_order {hole_hydro_(gz|dens|pfdens)}]
    if {[info exists state(hydro_scheme)] && $state(hydro_scheme) in {gz dens pfdens}} {
        catch {on_hydro_method_changed}
    }
    # Mean Profile coloured by a water scheme must ALSO refresh so it reflects the new
    # hydration (its stale files were just deleted, so this rebuilds fresh) - #12.
    # refresh_mean_surface_if_shown is a no-op unless the 3D surface is currently shown.
    if {[info exists state(mean_hydro_scheme)] && $state(mean_hydro_scheme) in {gz dens}} {
        catch {refresh_mean_surface_if_shown}
        catch {draw_mean_profile}
    }
    # Water schemes (G(z)/density) are now available - add them to the
    # property pickers, which otherwise hide them until hydration exists.
    catch {refresh_property_scheme_menus}
    set state(status) "Hydration: $total_w water-in-pore counts over $nfdata frame(s); selection \"$wsel\".$_anchor_note"
    vmdcon -info "VMDHole: hydration complete — $total_w water-in-pore counts over $nfdata frame(s) ($nfwater with water)."
    return 1
}

