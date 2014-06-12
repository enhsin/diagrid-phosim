package require Rappture

# open the XML file containing the run parameters
set driver [Rappture::library [lindex $argv 0]]

set installdir [file dirname [info script]]
set phosimDir [file join $installdir ..]
set exampleDir [file join $phosimDir examples]

set choice [$driver get input.(Mode).current]
if {$choice == "Example"} {
  set instanceCatalog [file join $exampleDir [$driver get input.group.(example).(instanceCatalog).current]]
  set extraCommands [file join $exampleDir [$driver get input.group.(example).(extraCommands).current]]
  set e2adc [$driver get input.group.(example).boolean(e2adc).current]
} elseif {$choice == "User"} {
  set instanceCatalogData [$driver get input.group.(user).(instanceCatalog).current]
  set instanceCatalog "instanceCatalog[pid]"
  set fid [open $instanceCatalog w]
  puts $fid $instanceCatalogData
  close $fid
  set extraCommandsData [$driver get input.group.(user).(extraCommands).current]
  set extraCommands "extraCommands[pid]"
  set fid [open $extraCommands w]
  puts $fid $extraCommandsData
  close $fid
  set e2adc [$driver get input.group.(user).boolean(e2adc).current]
} else {
  exit 0
} 
set binDir [file join $phosimDir bin]
set dataDir [file join $phosimDir data]
set sedDir [file join $phosimDir data SEDs]
set imageDir [file join $phosimDir data images]

file mkdir work output
set outputDir "output/"
set workDir "work/"

if {$e2adc == "yes"} {
    set e2adcflag 1
} else {
    set e2adcflag 0
}

set obsID "9999"
set pattern "obshistid"
set fid [open $instanceCatalog r]
while {[gets $fid line] >= 0} {
    if {[regexp $pattern $line]} {
       set obsID [lindex [split $line] end]
    }
}
close $fid

set status [catch {Rappture::exec python $binDir/phosim.py $instanceCatalog -c $extraCommands -e $e2adcflag -o $outputDir -w $workDir -b $binDir -d $dataDir --sed=$sedDir --image=$imageDir} out]

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
