# Tcl and headless use

VMDHole's pore workflow can be driven from the VMD Tk Console or from a script
passed to `vmd -dispdev text -e`. Pin the VMDHole version used by a script,
because state names and result fields can change between releases.

## Minimal pore script

```tcl
source /absolute/path/to/VMDHole/vmdhole/vmdhole.tcl
::VMDHole::init_executables

set molid [mol new /data/channel.pdb]
mol addfile /data/channel.xtc molid $molid waitfor all

set ::VMDHole::state(molid) $molid
set ::VMDHole::state(selection) "protein"
set ::VMDHole::state(frame_spec) "0:10:1000"
set ::VMDHole::state(cpoint) "12.3 4.5 -6.7"
set ::VMDHole::state(cvect) "0 0 1"
set ::VMDHole::state(radius_file) "/opt/hole2/rad/simple.rad"
set ::VMDHole::state(display_mode) "none"
set ::VMDHole::state(work_dir) "/data/results/channel"
set ::VMDHole::state(save_results) 1

if {![::VMDHole::run_analysis]} {
    puts stderr "VMDHole failed: $::VMDHole::state(status)"
    exit 1
}

foreach frame $::VMDHole::result_frames {
    set profile [dict get $::VMDHole::results $frame profile]
    puts "$frame,[dict get $profile min_radius]"
}
exit 0
```

Run it with:

```sh
vmd -dispdev text -e analyse.tcl
```

`run_analysis` returns `1` when the run completes and `0` for a cancelled or
failed run. Hard validation errors can also be raised as Tcl errors before
execution; production scripts should wrap the call in `catch` and return a
nonzero process status.

## Set native executable paths

`init_executables` reads saved configuration. Paths can also be set explicitly:

```tcl
set ::VMDHole::state(hole_exec) "/opt/vmdhole/bin/hole"
set ::VMDHole::state(sph_process_exec) "/opt/vmdhole/bin/sph_process"
set ::VMDHole::state(sos_triangle_exec) "/opt/vmdhole/bin/sos_triangle"
```

An empty path permits an embedded fallback where supported. This is useful for
portability but can be much slower.

## Frame specification

`state(frame_spec)` accepts the same syntax as the GUI: `now`, `all`, a frame
number, `start:end`, or `start:stride:end`.

## Read results without dialogs

GUI export procedures open file-selection dialogs and should not be called in
text mode. Read result dictionaries or metric helpers directly and write output
in the calling script:

```tcl
foreach frame $::VMDHole::result_frames {
    set metrics [::VMDHole::metrics_for_frame $frame]
    if {$metrics eq ""} { continue }
    puts "$frame,[dict get $metrics min_radius],[dict get $metrics volume]"
}
```

Record the VMDHole version with scripted output.

## Headless limitations

Do not call `vmdhole_tk`, `show_gui`, dialog procedures, or GUI export commands
under `-dispdev text`. Load coordinates before calling the analysis and ensure
that output paths are absolute and writable. Console logs contain engine paths,
warnings, failed frames, and the saved result root; capture them with the batch
job.

Tunnel analysis can also be configured through `state(tunnel_...)` and
`run_tunnel_analysis`, but the state keys are more extensive. Use the values in
the [parameter reference](parameters.md#tunnel-search), generate a run once in
the GUI, and pin that VMDHole version before automating a tunnel pipeline.
