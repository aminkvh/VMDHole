if {![package vsatisfies [package provide Tcl] 8.5]} {return}
package ifneeded vmdpathfinder 1.0.1 [list source [file join $dir vmdpathfinder.tcl]]
