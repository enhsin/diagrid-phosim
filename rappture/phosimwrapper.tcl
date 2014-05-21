package require Rappture

# open the XML file containing the run parameters
set driver [Rappture::library [lindex $argv 0]]

set instanceCatalog [$driver get input.(instanceCatalog).current]
set extraCommands [$driver get input.(extraCommands).current]
set e2adc [$driver get input.boolean(e2adc).current]
set outputDir "output/"
set obsID "99999999"

if {$e2adc == "yes"} {
    set e2adcflag 1
} else {
    set e2adcflag 0
}

set status [catch {Rappture::exec python phosim.py $instanceCatalog -c $extraCommands -e $e2adcflag} out]
$driver put output.log $out

if {$status == 0} {
    foreach fn [glob -nocomplain -types f -directory $outputDir *$obsID*] {
        set str [split $fn .]
        set fname [lindex [split [lindex $str {0}] /] end]
        set fileType "."
        append fileType [lindex $str {1}] "." [lindex $str {2}]

        set f [open $fn "r"]
        fconfigure $f -translation binary
        set data [read $f]

        $driver put output.string($fname).about.label $fname
        $driver put output.string($fname).current $data
        $driver put output.string($fname).filetype $fileType
        close $f
    }
}

# save the updated XML describing the run...
Rappture::result $driver
exit 0
