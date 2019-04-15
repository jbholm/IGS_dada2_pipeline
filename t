[1mdiff --git a/illumina_dada2.pl b/illumina_dada2.pl[m
[1mindex bf48cda..5836e61 100755[m
[1m--- a/illumina_dada2.pl[m
[1m+++ b/illumina_dada2.pl[m
[36m@@ -498,12 +498,6 @@[m [mif ( ( !@dbg ) || grep( /^qiime_and_validation$/, @dbg ) ) {[m
           find_raw_files( $inDir, $oneStep, $logFH );[m
     } else {[m
         if ($oneStep) {[m
[31m-[m
[31m-  # Correct??? Documentation say that r1 -> R1.fastq and r2, r3, r4 -> R2.fastq.[m
[31m-  # What actually happened was, R1 and R2 files autofound[m
[31m-  # from $inDir were used for barcodes, and silently nothing happened for[m
[31m-  # split_libraries (because the two files were already in the run directory)[m
[31m-  # pipeline was never coded to run 1-step WITHOUT -i[m
             $index1Input = $readsForInput = $r1file;[m
             $index2Input = $readsRevInput = $r2file;[m
         } else {[m
[36m@@ -524,61 +518,6 @@[m [mif ( ( !@dbg ) || grep( /^qiime_and_validation$/, @dbg ) ) {[m
 [m
     # In the multidbg branch, combine this section with qiime_and_validate[m
 [m
[31m-    my @cmds;[m
[31m-[m
[31m-    my %localNames;[m
[31m-    @localNames{ ( "readsFor", "readsRev", "index1", "index2" ) } =[m
[31m-      convert_to_local_if_gz( $wd, $readsForInput, $readsRevInput,[m
[31m-        $index1Input, $index2Input );[m
[31m-[m
[31m-    # Get forward and reverse reads as .fastq[m
[31m-    # if the index files aren't .gz, read in place. Otherwise...[m
[31m-[m
[31m-    if ( $readsForInput =~ /.gz$/ ) {[m
[31m-[m
[31m-        # overwrite any existing file in working dir[m
[31m-        # (we can't test if it's complete)[m
[31m-        my $dest = $localNames{"readsFor"};[m
[31m-        push( @cmds, "gzip --decompress --force < $readsForInput > $dest" )[m
[31m-          ;    #copy to our working dir[m
[31m-        $readsForInput = $dest;[m
[31m-        if ($oneStep) {[m
[31m-            $index1Input = $dest;[m
[31m-        }[m
[31m-    }[m
[31m-    if ( $readsRevInput !~ /.gz$/ ) {    # same for reverse reads[m
[31m-    } else {[m
[31m-        my $dest = $localNames{"readsRev"};[m
[31m-        push( @cmds, "gzip --decompress --force < $readsRevInput > $dest" );[m
[31m-        $readsRevInput = $dest;[m
[31m-        if ($oneStep) {[m
[31m-            $index2Input = $dest;[m
[31m-        }[m
[31m-    }[m
[31m-[m
[31m-    if ( !$oneStep ) {[m
[31m-[m
[31m-        # For two step, we need to do the same for the index files[m
[31m-        if ( $index1Input !~ /.gz$/ ) {[m
[31m-        } else {[m
[31m-            my $dest = $localNames{"index1"};[m
[31m-            push( @cmds, "gzip --decompress --force < $index1Input > $dest" );[m
[31m-            $index1Input = $dest;[m
[31m-        }[m
[31m-        if ( $index2Input !~ /.gz$/ ) {    # same for reverse reads[m
[31m-        } else {[m
[31m-            my $dest = $localNames{"index2"};[m
[31m-            push( @cmds, "gzip --decompress --force < $index2Input > $dest" );[m
[31m-            $index2Input = $dest;[m
[31m-        }[m
[31m-    }[m
[31m-[m
[31m-    # Execute the commands queued above[m
[31m-    if ( scalar @cmds > 0 ) {[m
[31m-        print "---Uncompressing raw file(s) to $wd\n";[m
[31m-        execute_and_log( @cmds, 0, $dryRun );[m
[31m-        @cmds = ();[m
[31m-    }[m
     if ($oneStep) {[m
         print $logFH[m
           "$project raw reads files:\n\t$readsForInput\n\t$readsRevInput\n";[m
[36m@@ -607,16 +546,50 @@[m [mif ( ( !@dbg ) || grep( /^extract_barcodes$/, @dbg ) ) {[m
 [m
     if ( !-e $barcodes ) {    # FIXME test if barcodes is the right size[m
         print "---Barcode files not found or were empty\n";[m
[31m-        my @cmds;[m
 [m
         my $start = time;[m
         print "---Extracting barcodes\n";[m
         $step1 = "extract_barcodes.py";[m
 [m
[31m-        my @rawfiles = find_raw_files( $inDir, $oneStep, $logFH );[m
[32m+[m[32m        my $index1Input[m
[32m+[m[32m          ;    # In one-step runs, these are the SAME files pointed to by[m
[32m+[m[32m        my $index2Input;    # $readsForInput and $readsRevInput[m
[32m+[m[32m        my $readsForInput;[m
[32m+[m[32m        my $readsRevInput;[m
[32m+[m[32m        ( $readsForInput, $readsRevInput, $index1Input, $index2Input ) =[m
[32m+[m[32m            $inDir ? find_raw_files( $inDir, $oneStep, $logFH )[m
[32m+[m[32m        : $oneStep ? ( $r1file, $r2file, $r1file, $r2file )[m
[32m+[m[32m        :            ( $r1file, $r4file, $r2file, $r3file );[m
         my %localNames;[m
[31m-        my ( $readsFor, $readsRev, $index1, $index2 ) =[m
[31m-          convert_to_local_if_gz( $wd, @rawfiles );[m
[32m+[m[32m        @localNames{ ( "readsFor", "readsRev", "index1", "index2" ) } =[m
[32m+[m[32m          convert_to_local_if_gz( $wd, $readsForInput, $readsRevInput,[m
[32m+[m[32m            $index1Input, $index2Input );[m
[32m+[m
[32m+[m[32m        # Get index files as .fastq[m
[32m+[m[32m        my @cmds;[m
[32m+[m[32m        if ( $index1Input !~ /.gz$/ ) {[m
[32m+[m
[32m+[m[32m            # if the index files aren't .gz, just read in place.[m
[32m+[m[32m        } else {    # Otherwise...[m
[32m+[m[32m            my $dest = $localNames{"index1"};[m
[32m+[m
[32m+[m[32m            # decompress to our project/run directory[m
[32m+[m[32m            push( @cmds, "gzip --decompress --force < $index1Input > $dest" );[m
[32m+[m[32m            $index1Input = $dest;[m
[32m+[m[32m        }[m
[32m+[m[32m        if ( $index2Input !~ /.gz$/ ) {    # same for reverse reads[m
[32m+[m[32m        } else {[m
[32m+[m[32m            my $dest = $localNames{"index2"};[m
[32m+[m[32m            push( @cmds, "gzip --decompress --force < $index2Input > $dest" );[m
[32m+[m[32m            $index2Input = $dest;[m
[32m+[m[32m        }[m
[32m+[m
[32m+[m[32m        # Execute the commands queued above[m
[32m+[m[32m        if ( scalar @cmds > 0 ) {[m
[32m+[m[32m            print "---Uncompressing raw index file(s) to $wd\n";[m
[32m+[m[32m            execute_and_log( @cmds, 0, $dryRun );[m
[32m+[m[32m            @cmds = ();[m
[32m+[m[32m        }[m
 [m
         my $mapOpt = "";[m
         my $bcLen  = 8;    # 2-step pcr[m
[36m@@ -628,7 +601,7 @@[m [mif ( ( !@dbg ) || grep( /^extract_barcodes$/, @dbg ) ) {[m
             $bcLen = 12;    # 1-step pcr[m
         }[m
         $cmd =[m
[31m-"extract_barcodes.py -f $index1 -r $index2 -c barcode_paired_end --bc1_len $bcLen --bc2_len $bcLen $mapOpt -o $wd";[m
[32m+[m[32m"extract_barcodes.py -f $index1Input -r $index2Input -c barcode_paired_end --bc1_len $bcLen --bc2_len $bcLen $mapOpt -o $wd";[m
         print "\tcmd=$cmd\n" if $verbose;[m
         print "---Waiting for barcode extraction to complete.\n";[m
 [m
[36m@@ -668,8 +641,17 @@[m [mif ( !@dbg || grep( /^demultiplex$/, @dbg ) ) {[m
     my @cmds;[m
     my $rForSeqsFq = "$rForSplit/seqs.fastq";[m
     my $rRevSeqsFq = "$rRevSplit/seqs.fastq";[m
[31m-    my ( $readsForInput, $readsRevInput, $index1Input, $index2Input ) =[m
[31m-      find_raw_files( $wd, $oneStep, $logFH );[m
[32m+[m
[32m+[m[32m    my $readsForInput;[m
[32m+[m[32m    my $readsRevInput;[m
[32m+[m[32m    my $index1Input;[m
[32m+[m[32m    my $index2Input;[m
[32m+[m[32m    ( $readsForInput, $readsRevInput, $index1Input, $index2Input ) =[m
[32m+[m[32m        $inDir ? find_raw_files( $inDir, $oneStep, $logFH )[m
[32m+[m[32m      : $oneStep ? ( $r1file, $r2file, $r1file, $r2file )[m
[32m+[m[32m      :            ( $r1file, $r4file, $r2file, $r3file );[m
[32m+[m
[32m+[m[32m# split_libraries_fastq.py accepts fastq.gz, so no need to convert_to_local_drop_gz[m
 [m
     # Just helps with printing[m
     my $revName;[m
[36m@@ -912,7 +894,6 @@[m [mif ( !@dbg || grep( /^demultiplex$/, @dbg ) ) {[m
     push( @cmds, "rm -rf $wd/${run}_R1.fastq" );[m
     push( @cmds, "rm -rf $wd/${run}_R2.fastq" );[m
     push( @cmds, "rm -rf $wd/${run}_R3.fastq" );[m
[31m-    push( @cmds, "rm -rf $wd/${run}_R2.fastq" );[m
     push( @cmds, "rm -rf $wd/${run}_I1.fastq" );[m
     push( @cmds, "rm -rf $wd/${run}_I2.fastq" );[m
     execute_and_log( @cmds, 0, $dryRun );[m
[36m@@ -1531,6 +1512,7 @@[m [msub find_raw_files {[m
     } else {[m
         my @i1s = glob("$wd/*I1.fastq $wd/*I1.fastq.gz");[m
         my @i2s = glob("$wd/*I2.fastq $wd/*I2.fastq.gz");[m
[32m+[m
         # also globbed for r1 and r2 files above[m
         my @r3s = glob("$wd/*R3.fastq $wd/*R3.fastq.gz");[m
         my @r4s = glob("$wd/*R4.fastq $wd/*R4.fastq.gz");[m
[36m@@ -1601,7 +1583,8 @@[m [msub convert_to_local_if_gz {[m
             my @dirs = File::Spec->splitdir($wd);[m
             $file =[m
                 "$wd/$dirs[scalar(@dirs) - 2]" . "_"[m
[31m-              . "$dirs[scalar(@dirs) - 1]" . "_$suffix";[m
[32m+[m[32m              . "$dirs[scalar(@dirs) - 1]"[m
[32m+[m[32m              . "_$suffix";[m
         }[m
         push( @ans, $file );[m
     }[m
[36m@@ -1777,6 +1760,7 @@[m [msub check_error_log {[m
 sub execute_and_log {[m
     my $dryRun = pop @_;[m
     my $logFH  = pop @_;[m
[32m+[m
     # CAREFUL, $cmd holds a reference to a variable in the caller![m
     foreach my $cmd (@_) {[m
         print "\t$cmd\n" if $verbose;[m
