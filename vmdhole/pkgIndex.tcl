if {![package vsatisfies [package provide Tcl] 8.5]} {return}
package ifneeded vmdhole 1.0.0 [list source [file join $dir vmdhole.tcl]]
