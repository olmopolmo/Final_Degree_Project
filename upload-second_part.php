<!DOCTYPE html>

<html>

<head>
    <!-- <link rel="stylesheet" type="text/css" href="styles.css"> -->
    <title>CRISPR-A Results</title>
    <meta http-equiv="content-type" content="text/html; charset=utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="description" content="" />
    <meta name="keywords" content="" />
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Public+Sans:wght@400;700&display=swap" rel="stylesheet">
    <script src="../../js/jquery.min.js"></script>
    <script src="../../js/slick.js"></script>
    <script src="../../js/scripts.js"></script>
    <script src="../../js/config.js"></script>
    <script src="../../js/skel.min.js"></script>
    <script src="https://demos.jquerymobile.com/1.4.2/js/jquery.js"></script>
	<script src='https://cdn.plot.ly/plotly-2.9.0.min.js'></script>
    <script src="https://cdn.jsdelivr.net/npm/igv@2.10.5/dist/igv.min.js"></script>
    <link rel="stylesheet" href="../../css/skel-noscript.css" />
    <link rel="stylesheet" href="../../css/normalise.css" />
    <link rel="stylesheet" href="../../css/slick.css" />
    <link rel="stylesheet" href="../../css/style.css" />
    <link rel="stylesheet" href="../../css/style-desktop.css" />
    <link rel="stylesheet" href="../../css/style-tablet.css" />
    <link rel="stylesheet" href="../../css/style-mobile.css" />
    <link rel="stylesheet" href="https://www.w3schools.com/w3css/4/w3.css">

</head>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.2/js/bootstrap.min.js"></script>
	<script>
        //Remove elements function
        function removeElement(arrayName1,arrayName2,arrayName3,arrayElement)
        {
            var len = arrayName1.length;
            for(var i=(len-1); i>=0 ;i-- )
            { 
                if(arrayName1[i]==arrayElement){
                    arrayName1.splice(i,1); 
                    arrayName2.splice(i,1);
                    arrayName3.splice(i,1);
                }
                    
            }
        };
		//GenomeBrowser function js
		function setReference(sample){
		var igvDiv = document.getElementById("igv-div");
		igv.removeAllBrowsers();
		var options = {

			reference: {
				id: sample,
				fastaURL:  "Output/Data/" + sample.concat("_reference-correctOrient.fasta"),
				indexURL:  "Output/Data/" + sample.concat("_reference-correctOrient.fasta.fai")
			},

			tracks: [
			{
				name: "test",
				type: "alignment",
				url: "Output/Data/" + sample  +".sorted.bam",
				indexURL: "Output/Data/" + sample + ".sorted.bam.bai",
				format: "bam",
				sort: {
					position: cut_site,
					option: "INSERT_SIZE",
					direction: "ASC"
				}
			}
			]
		};
		igv.createBrowser(igvDiv,options)
			.then(function(browser){
				console.log("Created IGV browser");
			})
		}
    function activateColapse(){
        // colapsible
            var coll = document.getElementsByClassName("collapsible");
            var i;
            for (i = 0; i < coll.length; i++) {
              coll[i].addEventListener("click", function() {
                this.classList.toggle("active");
                var content = this.nextElementSibling;
                if (content.style.display === "block") {
                  content.style.display = "none";
                } else {
                  content.style.display = "block";
                }
              });
            }
        }
    
	//function to create table and plot
	function PlotTable(num,metadels,metains, json_total_reads, cut_site,seq_len){
            var start = parseInt(document.getElementById("start_"+ num.toString()).value);
            var end = parseInt(document.getElementById("end_"+num.toString()).value);
            
           var remain=0;

            //new accumulative deletions
            var n = seq_len; //must be changed to the sequence length
            var del = Array(n).fill(0); //Array containing n 0 where n is the sequence length
            var count_dels = Array(n).fill(0); 
            var a = 1;
            var del_pos =  Array.from(Array(n), ()=>a++);
            for(l of metadels){
                if ((start<parseInt(l.Start) && parseInt(l.Start)<end) || (start<l.End && l.End<end)){
                    var s = parseInt(l.Start);
                    var e = l.End;
                    let change = Array.from(Array(e-s+1) , ()=>s++);
                    for (let v of change){
                        del[v] += ((parseInt(l.Freq))/(parseInt(l.changes))*100);
                        count_dels[v] += parseInt(l.Freq);
                    }
                    remain += l.Freq;
                }
            }

            //new accumulative deletions
            var n = seq_len; //must be changed to the sequence length
            var ins = Array(n).fill(0); //Array containing n 0 where n is the sequence length
            var count_ins = Array(n).fill(0); 
            var a = 1;
            var ins_pos =  Array.from(Array(n), ()=>a++);
            for(l of metains){
                if ((start<parseInt(l.Start) && parseInt(l.Start)<end) || (start<l.End && l.End<end)){
                    var s = parseInt(l.Start);
                    var e = l.End;
                    let change = Array.from(Array(e-s+1) , ()=>s++);
                    for (let v of change){
                        ins[v] += ((parseInt(l.Freq))/(parseInt(l.changes))*100);
                        count_ins[v] += parseInt(l.Freq);
                    }
                    remain += l.Freq;
                }
            }          
            removeElement(ins,ins_pos,count_ins,0)
            removeElement(del,del_pos,count_dels,0)


            var data_d = {
                x: del_pos,
                y: del,
                type: 'bar',
                hovertemplate: '<i>Percentage: </i>: %{y:.2f}' +
                '<br>Counts:<b>%{text}</b><extra></extra>',
                text: count_dels,
                showlegend: false,
                marker: {	
					color: 'rgb(147, 149, 247)'
				}
            };
        
            var data_i = {
                x: ins_pos,
                y: ins,
                xaxis: 'x2',
                yaxis: 'y2',
                type:'bar',
                text: count_ins,
                hovertemplate: '<i>Percentage: </i>: %{y:.2f}' +
                '<br>Counts:<b>%{text}</b><extra></extra>',
                showlegend: false,
                marker: {	
					color: 'rgb(147, 149, 247)'
				}
            };    
        
            var data = [data_d,data_i];
            if(start<cut_site && cut_site<end){
                var layout = {
                showlegend: false,
                shapes: [{
                    type: "line",
                    y0:0,
                    y1:100,
                    x0:cut_site,
                    x1:cut_site,
                    yref: "paper",
                    line: {color:"red",dash:"dot"},
                    name: "cut site"
                    },
                    {
                    type: "line",
                    y0:0,
                    y1:100,
                    x0:cut_site,
                    x1:cut_site,
                    yref: "paper",
                    xref: "x2",
                    line: {color:"red",dash:"dot"},
                    name: "cut site"
                    }
                ],
                annotations: [
                    {
                    x: cut_site,
                    xref: 'X',
                    yref: 'paper',
                    text: 'cut site'
                    },
                    {
                    x: cut_site,
                    xref: 'x2',
                    yref: 'paper',
                    text: 'cut site'
                    }
                ],
                grid: {rows: 1, columns: 2, pattern: 'independent'},
                yaxis: {
                    title: {
                    text: 'Reads Percentage(%)',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                },
                xaxis: {
                    title: {
                    text: 'Accumulative Deletions Locations',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                },
                xaxis2: {
                    title: {
                    text: 'Accumulative Insertions Locations',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                },
                yaxis2: {
                    title: {
                    text: 'Reads Percentage(%)',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                }
                };

            } else{
                var layout = {
                showlegend: false,
                grid: {rows: 1, columns: 2, pattern: 'independent'},
                yaxis: {
                    title: {
                    text: 'Reads Percentage(%)',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                },
                xaxis: {
                    title: {
                    text: 'Accumulative Deletions Locations',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                },
                xaxis2: {
                    title: {
                    text: 'Accumulative Insertions Locations',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                },
                yaxis2: {
                    title: {
                    text: 'Reads Percentage(%)',
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                    }
                }
                };
            }
            Plotly.newPlot('myplo_'+num.toString(), data, layout);
            var values = [del_pos,del,ins_pos,ins]

            var data = [{
            type: 'table',
            header: {
                values: [["<b>Position Deletions</b>"],["<b>Percentage Deletions</b>"],
                ["<b>Position Insertions</b>"],["<b>Percentage Insertions</b>"]],
                align: "center",
                line: {width: 1, color: 'black'},
                fill: {color: "grey"},
                font: {family: "Arial", size: 12, color: "white"}
            },
            cells: {
                values: values,
                align: "center",
                line: {color: "black", width: 1},
                font: {family: "Arial", size: 11, color: ["black"]}
            }
            }];
            var layout = {}

            Plotly.newPlot('mytab_'+num.toString(), data);

            //read json for the moment I set the value
            var per = document.getElementById("percent_"+num.toString());
			var message = "<center>Modification Percentage = "
            per.innerHTML = message.concat(Math.round((remain/json_total_reads)*100*100)/100,"%","</center>")

        }
    </script>
    
    <style media="screen">	
        /**	
     * Sticky navigation	
     */	
    .sticky {	
    	background-color: #ffffff;	
    	position: -webkit-sticky;	
    	position: sticky;	
    	top: 0;	
        z-index: 600;	
    }	
    .collapsible {	
      background-color: white;	
      color: black;	
      cursor: pointer;	
      padding: 6 px;	
      width: 100%;	
      border: none;	
      text-align: left;	
      /* outline: none; */	
      font-size: 15px;	
    }	
    .active, .collapsible:hover {	
      background-color: #727272;	
    }	
    .content {	
      /* padding: 0 18px; */	
      display: none;	
      overflow: hidden;	
      background-color: #f1f1f1;	
    }	
    </style>	


<body class="results">

    <header class="header">
        <div class="header-logo">
            <a href="https://synbio.upf.edu/crispr-a/"><img src="../../Images/crispr-a-logo.svg" alt="CRISPR-A"></a>
        </div>
        <div class="header-nav">
            <a class="header-nav__link button button--black" href="../../documentation.html">Help</a>
            <a class="header-nav__link button button--white" href="">Reference</a>
        </div>
    </header>

    <!-- Progress bar holder -->
    <div id='progress' style='width:500px;border:1px solid #ccc; display: none;'></div>
    <!-- Progress information -->
    <div id='information'></div>

    <?php
    error_reporting(E_ALL);
    ini_set("display_errors", 1);
    //Include functions!
    include '../../functions.php';

    // retrieving samples array
    $file = 'samples.php';
    eval('$allSamples = ' . file_get_contents($file) . ';');

	// Run Nextflow
        //$cmd='nextflow run ../../crispr-ga_nextflow/crispr-ga_umis_subs.nf -c nextflow_input.config 2>&1';
        //$cd = chdir($target_dir);

	    // if ($cd) {
	    //    $cmdOutput = exec($cmd, $output, $var);
	    //    echo '<pre>'; print_r($output); echo '</pre>';
	    // }

        // If mocks.php exists the sample naming gets the preference
        if (file_exists("mocks.php")){
            $file_m = 'mocks.php';
            eval('$allmocks = ' . file_get_contents($file_m) . ';');
            }
    	//Number of samples
    	$countfiles1 = count($allSamples);

	    //Menu bar for sample results
        echo "	
        <nav class='sticky'>";
	    echo "
	    <nav class='samples-nav'>
	        <span class='samples-nav__label text--bold'>Shortcuts </span>";
	        // Check if there are less than 4 samples
	        if($countfiles1 < 4){
	            $shortMenu = $countfiles1;
	        } else {
	            $shortMenu = 4;
	        }
	        //Menu row with 4 or less samples
	        menuRow($allSamples, 0, $shortMenu, "samples-nav__list");

	        // Menu part for more than 4 samples
	        if($countfiles1 > 4){
	            echo "
	            <span class='samples-nav__button'>
	                <img src='../../Images/arrow-down.svg' alt=''>
	            </span>";
	            //menuRow($allSamples, 4, 8, "samples-nav__list--secondary");
	            menuRow($allSamples, 4, $countfiles1, "samples-nav__list--secondary");
	        }

	        echo "
	        </nav>";
        echo "	
        </nav>";

	    //Show results
	    echo "<main><div class='main__inner'>
	    <h1 class='headline'>CRISPR-A<br/>Results</h1>";

	    //Table with read counts summary after each data pre-processing step
	    echo "<h2 class='results-content__heading heading'>Summary</h2><div>";	
        echo "<h3 class='table__heading'>Pre-processing (reads count)</h3>";	
        echo "<div>";	
        echo "

	    <table class='results-table'>
	        <tr>
	            <th class='text text--bold'>Sample</th>
	            <th class='text text--bold'>Raw</th>
	            <th class='text text--bold'>Merged</th>
                <th class='text text--bold'>Adapters</th>
	            <th class='text text--bold'>Quality filtered</th>
	            <th class='text text--bold'>Clustered</th>
	            <th class='text text--bold'>Aligned</th>
	        </tr>";

	    for($i=0;$i<$countfiles1;$i++){
            if (file_exists("mocks.php")) {
                $filename="Output/Data/{$allmocks[$i]}_summary.csv";
            } else {
                $filename="Output/Data/sample_{$i}_summary.csv";
            }
	        $sample_name=preg_split("/[_]+/", $allSamples[$i]);
	        if (file_exists($filename)){
				echo "<tr> \n";
	             if (file_exists("mocks.php")){
                        $sample = $allmocks[$i]; //genome browser
                        //read cut site
                        $cut = file_get_contents("Output/Results/{$sample}_cutSite.json");
                        $cut_site = json_decode($cut,true);
                        if (empty($cut_site)){
                            echo"<script>var cut_site = 0;</script>"; //pass the variable to the javascript enviroment
                        }else{
                            echo"<script> var cut_site = $cut_site[0]; </script>";
                        }
                        }
                    else{
                         $sample = 'sample_' . $i; //genome browser
                        //read cut site
                        $cut = file_get_contents("Output/Results/{$sample}_cutSite.json");
                        $cut_site = json_decode($cut,true);
                        if (empty($cut_site)){
                            echo"<script>var cut_site = 0;</script>"; //pass the variable to the javascript enviroment
                        }else{
                            echo"<script> var cut_site = $cut_site[0]; </script>";
                        }
                        }
                echo"<script> var sample_$i = '$sample'</script>"; //genome browser
                echo"<td><button class='cta button button--black' ref='#' onClick='setReference(sample_$i)'>$sample_name[0]</button></td>"; //genome browser
				//echo "<td class='text'> $sample_name[0] </td>";
	            $f_summary=fopen($filename,"r"); // file with read counts
	            $header=0;
	            while (($line = fgetcsv($f_summary)) !== false) {
	                if ( $header == 0 ){
	                    $header = 1;
	                } else {
	                    echo "          <td class='text'>" . $line[1] . " </td> \n";
	                }
	            }
	            echo "</tr> \n";
	            fclose($f_summary);
	        }
	    }
    	//start of genome browser php part
        echo '<div id="igv-div" style="padding-top: 10px;padding-bottom: 10px; border:1px solid lightgray"></div>';

	    echo "
	    </table
        </div>
	    </div>";
        echo "<br>
        <br>";
	    //Table with edits summary
        echo "<button type='button' class='collapsible'><h3 class='table__heading'>Edition classes (%)</h3></button>";	
        echo "<div class='content'>";
	    echo "<div class='table__wrapper'>

	    <table class='results-table'>
	        <tr>
	            <th class='text text--bold'>Sample</th>
	            <th class='text text--bold'>Wt</th>
	            <th class='text text--bold'>Template-based</th>
	            <th class='text text--bold'>Indels</th>
                <th class='text text--bold'>Delins</th>
	            <th class='text text--bold'>Insertions</th>
	            <th class='text text--bold'>Ins inframe</th>
	            <th class='text text--bold'>Ins outfarme</th>
	            <th class='text text--bold'>Deletions</th>
	            <th class='text text--bold'>Dels inframe</th>
	            <th class='text text--bold'>Dels outfarme</th>
	        </tr>";

	    for($i=0;$i<$countfiles1;$i++){
            if (file_exists("mocks.php")) {
                $filename="Output/Results/{$allmocks[$i]}_edits.csv";
            } else {
                $filename="Output/Results/sample_{$i}_edits.csv";
            }
	        if (file_exists($filename)){
	            echo "<tr> \n";
	            $sample_name=preg_split("/[_]+/", $allSamples[$i]);
	            echo "<td class='text'> $sample_name[0] </td>";
	            $f_summary=fopen($filename,"r"); // file with edit counts
	            $header=0;
	            while (($line = fgetcsv($f_summary)) !== false) {
	                if ( $header == 0 ){
	                    $header = 1;
	                } else {
	                    echo "          <td class='text'>" . $line[2] . " </td> \n";
	                }
	            }
	            echo "</tr> \n";
	            fclose($f_summary);
	        }
	    }

        echo "
        </table
        </div>
        </div>
        </div>";

        if (file_exists("mocks.php")){
            //Table with editions after mock filtering
            echo "<button type='button' class='collapsible'><h3 class='table__heading'>Noise Substraction</h3></button>";	
            echo "<div class='content'>";
    	    echo "<div class='table__wrapper'>

    	    <table class='results-table'>
    	        <tr>
    	            <th class='text text--bold'>Sample</th>
    	            <th class='text text--bold'>Denoised reads</th>
                    <th class='text text--bold'>Indels Denoised</th>
                    <th class='text text--bold'>Indels In Peak</th>
                    <th class='text text--bold'>Indels Out Peak</th>
                    <th class='text text--bold'>Deletions Denoised</th>
                    <th class='text text--bold'>Deletions In Peak</th>
                    <th class='text text--bold'>Deletions Out Peak</th>
    	            <th class='text text--bold'>Noisy Reads</th>
    	            <th class='text text--bold'>Noisy Positions</th>
                    <th class='text text--bold'>Normalization Ratio</th>
    	        </tr>";

            for ($i = 0; $i < count($allmocks); $i++) {
                $searchword = 'sample_' . $i;
                $matches = array_filter($allmocks, function($var) use ($searchword) { return preg_match("/\b$searchword\b/i", $var); });
                $indexes = (array_keys($matches));
                if (count($matches) != 0){
                    for ($j = 0; $j < count($indexes); $j++) {
                        if ($allmocks[$indexes[$j]] != $allmocks[$indexes[0]]) {
                	        $filename="Output/Results/{$allmocks[$indexes[$j]]}_MockSummary.csv";
                	        if (file_exists($filename)){
                	            echo "<tr> \n";
                                $sample_name=preg_split("/[_]+/", $allSamples[$indexes[$j]]);
                	            echo "<td class='text'> $sample_name[0] </td>";
                	            $f_summary=fopen($filename,"r"); // file with edit counts
                	            $header=0;
                	            while (($line = fgetcsv($f_summary)) !== false) {
                	                if ( $header == 0 ){
                	                    $header = 1;
                	                } else {
                	                    echo "          <td class='text'>" . $line[2] . " </td> \n";
                	                }
                	            }
                                echo "</tr> \n";
                                fclose($f_summary);
                            }
            	        }
            	    }
                }
            }
            echo "
            </table
            </div>
            </div>
            </div>";

            //Table with edits summary
            echo "<button type='button' class='collapsible'><h3 class='table__heading'>Edition Classes After Denoise (%)</h3></button>";	
            echo "<div class='content'>";
    	    echo "<div class='table__wrapper'>

            <table class='results-table'>
    	        <tr>
    	            <th class='text text--bold'>Sample</th>
    	            <th class='text text--bold'>Wt</th>
    	            <th class='text text--bold'>Template-based</th>
    	            <th class='text text--bold'>Indels</th>
    	            <th class='text text--bold'>Insertions</th>
    	            <th class='text text--bold'>Ins inframe</th>
    	            <th class='text text--bold'>Ins outfarme</th>
    	            <th class='text text--bold'>Deletions</th>
    	            <th class='text text--bold'>Dels inframe</th>
    	            <th class='text text--bold'>Dels outfarme</th>
    	        </tr>";

            for ($i = 0; $i < count($allmocks); $i++) {
                $searchword = 'sample_' . $i;
                $matches = array_filter($allmocks, function($var) use ($searchword) { return preg_match("/\b$searchword\b/i", $var); });
                $indexes = (array_keys($matches));
                if (count($matches) != 0){
                    for ($j = 0; $j < count($indexes); $j++) {
                        if ($allmocks[$indexes[$j]] != $allmocks[$indexes[0]]) {
                	        $filename="Output/Results/{$allmocks[$indexes[$j]]}_Mockedits.csv";
                	        if (file_exists($filename)){
                	            echo "<tr> \n";
                                $sample_name=preg_split("/[_]+/", $allSamples[$indexes[$j]]);
                	            echo "<td class='text'> $sample_name[0] </td>";
                	            $f_summary=fopen($filename,"r"); // file with edit counts
                	            $header=0;
                	            while (($line = fgetcsv($f_summary)) !== false) {
                	                if ( $header == 0 ){
                	                    $header = 1;
                	                } else {
                	                    echo "          <td class='text'>" . $line[2] . " </td> \n";
                	                }
                	            }
                	            echo "</tr> \n";
                	            fclose($f_summary);
                            }
                        }
                    }
    	        }
    	    }
           echo "
            </table
            </div>
            </div>";
        }
        echo "     
        </div>";
        
        echo "<script>activateColapse()</script>";
	    //Heatmap of indel length and position distances between samples when there are more than 2 samples
		if ($countfiles1 > 2) {
	        echo "<div class='sample__wrapper'>";
	        echo "<h2 class='heading'>Comparison - samples distance - </h2>";
	        echo "<h3 class='text'>Indels position and length Jensen distance</h3>";
            echo "<embed src='Output/Results/jensen_pos-len.html' width='1600' height='750'>";
	        echo "</div>";
	    }

        if (file_exists("mocks.php")){
            for ($i = 0; $i < count($allmocks); $i++) {
                $searchword = 'sample_' . $i;
                $matches = array_filter($allmocks, function($var) use ($searchword) { return preg_match("/\b$searchword\b/i", $var); });
                $indexes = (array_keys($matches));
                if (count($matches) != 0){
                    for ($j = 0; $j < count($indexes); $j++) {
                        $sample_name=preg_split("/[_]+/", $allSamples[$indexes[$j]]);

                        echo " <div class='sample__wrapper' id='$sample_name[0]'><h2 class='headline'>Sample $sample_name[0]</h2>";

                        echo " <h3 class='heading'>Analysis summary</h3>";

                        echo "<div class='results__blocks'>";
                        echo "<div class='results__block--third'>";
                        echo "<h4 class='text'>Pre-processing</h4>";
                        echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_reads.html' width='500' height='500'> ";
                        echo "</div>";
                        echo "<div class='results__block--third'>";
                        echo " <h4 class='text'>Indels quality</h4>";
                        echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_QC-indels.html' width='500' height='500'> ";
                        echo "</div>";
                        echo "<div class='results__block--third'>";
                        echo " <h4 class='text'>Edition</h4>";
                        echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_edition.html' width='500' height='500'> ";
                        echo "</div>";
                        echo "</div>
                              <br><br><br><br>";

                        if ($allmocks[$indexes[$j]] != $allmocks[$indexes[0]]) {
                            echo " <h3 class='heading'>Noise Substraction</h3>";

                            echo "<div class='results__blocks results__blocks--third'>";
                            echo "<div class='w3-container w3-half'>";
                            echo "<h4 class='text'>Treatment vs Mock deletions</h4>";
                            echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_MvsTdel.html' width='770' height='750'> ";
                            echo "</div>";
                            echo "<div class='w3-container w3-half'>";
                            echo " <h4 class='text'>Deletions corrected</h4>";
                            echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_delsSubstraction.html' width='770' height='750'> ";
                            echo "</div>";
                            echo "</div>";
                            echo "<div class='results__blocks results__blocks--third'>";
                            echo "<div class='w3-container w3-half'>";
                            echo " <h4 class='text'>Corrected Plot</h4>";
                            echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_MvsTins.html' width='770' height='750'> ";
                            echo "</div>";
                            echo "<div class='w3-container w3-half'>";
                            echo " <h4 class='text'>Insertions corrected</h4>";
                            echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_insSubstraction.html' width='770' height='750'> ";
                            echo "</div>";
                            echo "</div>
                            <br><br><br>";

                        }

                        echo " <h3 class='heading'>Editing results</h3>";
						$k = $i*100 + $j; //this variable is to keep track of the plots for accumulative

                        //read cut site
                        $cut = file_get_contents("Output/Results/{$allmocks[$indexes[$j]]}_cutSite.json");
                        $cut_site = json_decode($cut,true);
                        if (empty($cut_site)){
                            echo"<script>var cut_site_$k = 0;</script>"; //pass the variable to the javascript enviroment
                        }else{
                            echo"<script> var cut_site_$k = $cut_site[0]; </script>";
                        }

                        //read sequece length
                        $seq = file_get_contents("Output/Results/{$allmocks[$indexes[$j]]}_length.json");
                        $seq_len = json_decode($seq,true);
                        if (empty($seq_len)){
                            echo"<script>var seq_len_$k = 0;</script>"; //pass the variable to the javascript enviroment
                            $seq_len[0] = 0;
                        }else{
                            echo"<script> var seq_len_$k = $seq_len[0]; </script>";
                        }

                        echo " <h4 class='text'>Deletions</h4>";
                        echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_Deletions.html' width='1600' height='700'> ";
						echo " <h4 class='text'>Insertions</h4>";
                        echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_Insertions.html' width='1600' height='700'> ";
						echo " <h4 class='text'>Accumulative</h4>
    						<div class='form__block form__block--ngs'>
							<div class='form__wrapper' id='row1'>

								<div class='input__wrapper'>
									<br>
									<div class='input__inner input__inner--acc'>
										<label for='range-1a' class='form__label form__label--full'>Start:</label>
										<input class='input__text required' type='number' name='start' id='start_$k' min='0' max='10000000' value=0 /><br>
									</div>

									<div class='input__inner input__inner--acc'>
										<label  for='range-1b' class='form__label form__label--full'>End:</label>
										<input class='input__text'type='number' name='end' id='end_$k' min='0' max='10000000' value='$seq_len[0]'><br>
										<br>
									</div>
									&nbsp;

									<div class ='container'>
									<div class='center'>
										<button class='cta button button--black' onclick='PlotTable($k, metadel_$k,metains_$k, json_total_reads_$k, cut_site_$k,seq_len_$k)'>Refresh</button>	   
									</div>
									</div>

								</div>

							</div>
							</div>

						        <div id='percent_$k'  class='percent'></div>
                                <br>
                                <br>
                               <p id='myplo_$k'></p>
							   <p id='mytab_$k'></p>	
										   
						";
                        
                        //deletion  metadata json
                        if(file_exists("Output/Results/{$allmocks[$indexes[$j]]}_metadels.json")){
                            $data_d = file_get_contents("Output/Results/{$allmocks[$indexes[$j]]}_metadels.json");
                            echo"<script>
                            var metadel_$k = $data_d
                            </script>";
                        }else{
                            echo"<script>
                            var metadel_$k = []
                            </script>";
                        };
  
						//insertion metadata json
                        if(file_exists("Output/Results/{$allmocks[$indexes[$j]]}_metains.json")){
                            $data_i = file_get_contents("Output/Results/{$allmocks[$indexes[$j]]}_metains.json");
                            echo"<script>
                            var metains_$k =  $data_i
                            </script>";
                        }else{
                            echo"<script>
                            var metains_$k = []
                            </script>";
                        };

						//total reads json
                        if(file_exists("Output/Results/{$allmocks[$indexes[$j]]}_Total_reads.json")){
                            $total_reads = file_get_contents("Output/Results/{$allmocks[$indexes[$j]]}_Total_reads.json");
                            echo"<script>
                            var json_total_reads_$k = $total_reads
                            </script>";
                        }else{
                            echo"<script>
                            var json_total_reads_$k = []
                            </script>";
                        };

                        echo"	
                        <p id='error_acc_$k'></p>
                        <script>if(metadel_$k != [] || metains_$k != [] || json_total_reads_$k != [] || cut_site_$k != [] || seq_len_$k!=[]){
                                    PlotTable($k, metadel_$k,metains_$k, json_total_reads_$k,cut_site_$i,seq_len_$k)    
                                } 
                                else{
                                    document.getElementById('error_acc_$k').innerHTML ='No accumulative plot generated';
                                }
                        </script>";  
                        
                        echo " <h4 class='text'>Top alleles</h4>";
                        echo "<div class='results__blocks results__blocks--third'>";
                        echo "<div class='results__block--half'>";
                        echo " <embed type='text/html' src='Output/Results/{$allmocks[$indexes[$j]]}_top.html' width='750' height='750'> ";
                        echo "</div>";
                        echo "<div class='results__blocks--half'>";
                        echo " <div class='results__blocks results__blocks--full results-content__image'><img src='Output/Results/{$allmocks[$indexes[$j]]}_top-alleles_LOGO.png' width='750' height='750'></div>";
                        echo "</div>";
                        echo "</div>";

                        echo " <h4 class='text'>Substitutions</h4>";
                        echo " <div class='results__blocks results__blocks--full results-content__image'><img src='Output/Results/{$allmocks[$indexes[$j]]}_subs-perc_plot.png' width='100%'></div>";
						echo " <div class='results__blocks results__blocks--full results-content__image'><img src='Output/Results/{$allmocks[$indexes[$j]]}_subs-perc_plot_LOGO.png' width='100%'></div>";


                        echo "<h3 class='heading'>Edits complete information</h3>";
                        echo '<div class="input__wrapper">';
                        echo '<div class="input__inner input__inner--third">';	
                        echo "<span style='text-align: center;' class='form__label'>$sample_name[0]_subs-percentages</span>";
                        echo " <br> ";
                        echo "<a href='Output/Results/{$allmocks[$indexes[$j]]}_subs-perc.csv' download='$sample_name[0]_subs-perc.csv' class='button button--black download'>Download</a>";
                        echo " <br> ";
                      	                        echo '</div>';	
                        echo '<div class="input__inner input__inner--third">';	
                        echo "<span style='text-align: center;' class='form__label'>$sample_name[0]_indels</span>";	
                        echo " <br> ";	
                        echo "<a href='Output/Results/{$allmocks[$indexes[$j]]}_indels.csv' download='$sample_name[0]_indels.csv' class='button button--black download'>Download</a>";	
                        echo '</div>';	
                        echo '<div class="input__inner input__inner--third">';	
                        echo "<span style='text-align: center;' class='form__label'>$sample_name[0]_metadata</span>";	
                        echo " <br> ";	
                        echo "<a href='Output/Results/{$allmocks[$indexes[$j]]}_metadata.csv' download='$sample_name[0]_metadata.csv' class='button button--black download'>Download</a>";	
                        echo '</div>';

                        if ($allmocks[$indexes[$j]] != $allmocks[$indexes[0]]) {	
                            echo '<div class="input__inner input__inner--third">';	
                            echo "<span style='text-align: center;' class='form__label'>$sample_name[0]_indels_Correctedindels</span>";	
                            echo " <br> ";	
                            echo "<a href='Output/Results/{$allmocks[$indexes[$j]]}_CorrectedIndels.csv' download='$sample_name[0]_CorrectedIndels.csv' class='button button--black download'>Download</a>";	
                            echo '</div>';	
                        }	
                        echo '</div>';	
                        echo '</div>';	
                        }	
                    }	
                }	
            }
            else {
                for($i=0;$i<$countfiles1;$i++){

                    $sample_name=preg_split("/[_]+/", $allSamples[$i]);

                    if (file_exists("Output/Results/sample_{$i}_cutSite.json")){
                   	    //read cut site
                        $cut = file_get_contents("Output/Results/sample_{$i}_cutSite.json");
                        $cut_site = json_decode($cut,true);
    
                        if (empty($cut_site)){
                            echo"<script>var cut_site_$i = 0;</script>"; //pass the variable to the javascript enviroment
                        }else{
                            echo"<script> var cut_site_$i = $cut_site[0]; </script>";
                        }

                    }else{
                         echo"<script>var cut_site_$i = 0;</script>"; //pass the variable to the javascript enviroment
                    }
                    //read sequece length
                    $seq = file_get_contents("Output/Results/sample_{$i}_length.json");
                    $seq_len = json_decode($seq,true);
                    if (empty($seq_len)){
                        echo"<script>var seq_len_$i = 0;</script>"; //pass the variable to the javascript enviroment
                        $seq_len[0] = 0;
                    }else{
                        echo"<script> var seq_len_$i = $seq_len[0]; </script>";
                    }
                    
                    echo " <div class='sample__wrapper' id='$sample_name[0]'><h2 class='headline'>Sample $sample_name[0]</h2>";

                    echo " <h3 class='heading'>Analysis summary</h3>";

                    echo "<div class='results__blocks'>";
                    echo "<div class='results__block--third'>";
                    echo "<h4 class='text'>Pre-processing</h4>";
                    echo " <embed type='text/html' src='Output/Results/sample_{$i}_reads.html' width='500' height='500'> ";
                    echo "</div>";
                    echo "<div class='results__block--third'>";
                    echo " <h4 class='text'>Indels quality</h4>";
                    echo " <embed type='text/html' src='Output/Results/sample_{$i}_QC-indels.html' width='500' height='500'> ";
                    echo "</div>";
                    echo "<div class='results__block--third'>";
                    echo " <h4 class='text'>Edition</h4>";
                    echo " <embed type='text/html' src='Output/Results/sample_{$i}_edition.html' width='500' height='500'> ";
                    echo "</div>";
                    echo "</div>
                          <br><br><br><br>";

                    echo " <h3 class='heading'>Editing results</h3>";

                    echo " <h4 class='text'>Deletions</h4>";
                    echo " <embed type='text/html' src='Output/Results/sample_{$i}_Deletions.html' width='1600' height='700'> ";
					echo " <h4 class='text'>Insertions</h4>";
                    echo " <embed type='text/html' src='Output/Results/sample_{$i}_Insertions.html' width='1600' height='700'> ";
					echo " <h4 class='text'>Accumulative</h4>
   						<div class='form__block form__block--ngs'>
						<div class='form__wrapper' id='row1'>

							<div class='input__wrapper'>
								<br>
								<div class='input__inner input__inner--acc'>
									<label for='range-1a' class='form__label form__label--full'>Start:</label>
									<input class='input__text required' type='number' name='start' id='start_$i' min='0' max='10000000' value=0 /><br>
								</div>

								<div class='input__inner input__inner--acc'>
									<label  for='range-1b' class='form__label form__label--full'>End:</label>
									<input class='input__text'type='number' name='end' id='end_$i' min='0' max='10000000' value='$seq_len[0]'><br>
									<br>
								</div>
								&nbsp;
								&nbsp;
								<div class ='container'>
								<div class='center'>
									<button class='cta button button--black' onclick='PlotTable($i, metadel_$i,metains_$i, json_total_reads_$i,cut_site_$i,seq_len_$i)'>Refresh</button>	   
								</div>
								</div>

							</div>

						</div>
						</div>

                            <div id='percent_$i'  class='percent'></div>
                            <br>
                            <br>
                            <p id='myplo_$i'></p>
							<p id='mytab_$i'></p>

					";
                    
                //deletion metadata json
                if(file_exists("Output/Results/sample_{$i}_metadels.json")){
                    $data_d = file_get_contents("Output/Results/sample_{$i}_metadels.json");
                    echo"<script>
                    var metadel_$i = $data_d
                    </script>";
                }else{
                    echo"<script>
                    var metadel_$i = []
                    </script>";
                };
            
                //insertion metadata json
                if(file_exists("Output/Results/sample_{$i}_metains.json")){
                    $data_i = file_get_contents("Output/Results/sample_{$i}_metains.json");
                    echo"<script>
                    var metains_$i =  $data_i
                    </script>";
                }else{
                    echo"<script>
                    var metains_$i = []
                    </script>";
                };  
                //total reads json
                if(file_exists("Output/Results/sample_{$i}_Total_reads.json")){
                    $total_reads = file_get_contents("Output/Results/sample_{$i}_Total_reads.json");
                    echo"<script>
                    var json_total_reads_$i = $total_reads
                    </script>";
                }else{
                    echo"<script>
                    var json_total_reads_$i = []
                    </script>";
                };

                echo"	
                <p id='error_acc_$i'></p>
                <script>if(metadel_$i != [] || metains_$i != [] || json_total_reads_$i != [] || cut_site_$i != [] || seq_len_$i!=[]){
                            PlotTable($i, metadel_$i,metains_$i, json_total_reads_$i,cut_site_$i,seq_len_$i)    
                        } 
                        else{
                            document.getElementById('error_acc_$i').innerHTML ='No accumulative plot generated';
                        }
                </script>";  
					echo " <h4 class='text'>Top alleles</h4>";
                    echo "<div class='results__blocks results__blocks--third'>";
                    echo "<div class='results__block--half'>";
                    echo " <embed type='text/html' src='Output/Results/sample_{$i}_top.html' width='750' height='750'> ";
                    echo "</div>";
                    echo "<div class='results__blocks--half'>";
                    echo " <div class='results__blocks results__blocks--full results-content__image'><img src='Output/Results/sample_{$i}_top-alleles_LOGO.png' width='750' height='750'></div>";
                    echo "</div>";
                    echo "</div>";
                    echo " <h4 class='text'>Substitutions</h4>";
                    echo " <div class='results__blocks results__blocks--full results-content__image'><img src='Output/Results/sample_{$i}_subs-perc_plot.png' width='100%'></div>";
                    echo " <div class='results__blocks results__blocks--full results-content__image'><img src='Output/Results/sample_{$i}_subs-perc_plot_LOGO.png' width='100%'></div>";                    
                    echo "<h3 class='heading'>Edits complete information</h3>";
                    echo '<div class="input__wrapper">';
                    echo '<div class="input__inner input__inner--half">';
                    echo " <br> ";
                    echo "<span class='form__label'>$sample_name[0]_subs-percentages</span>";
                    echo "<a href='Output/Results/sample_{$i}_subs-perc.csv' download='$sample_name[0]_subs-perc.csv' class='button button--black download'>Download</a>";
                    echo " <br> ";
                    echo " <br> ";
                    echo "<span class='form__label'>$sample_name[0]_indels</span>";
                    echo "<a href='Output/Results/sample_{$i}_indels.csv' download='$sample_name[0]_indels.csv' class='button button--black download'>Download</a>";
                    echo '</div>';
                    echo '</div>';
                    echo '</div>';
                }
           }
	    ?>
	    </div>

	</main>

	<section class="footer-hero hero">
	    <div class="hero__text">
	        <h4 class="footer-hero__headline headline type--white">Reference </h4>
	        <h5 class="footer-hero__subheadline subheadline type--white">Guell M, Yang L, Church G (2014) </h5>
	        <h5 class="footer-hero__subheadline subheadline type--white">Genome Editing Assessment using CRISPR Genome Analyzer (CRISPR-GA) Bioinformatics.</h5>
	            <a class="footer-hero__button button button--green" href="">Read</a>
	    </div>
	</section>

	<footer class="footer">

	    <div class="footer__row">
	        <div class="footer__column">
	            <h4 class="footer__heading type--purple">Data controller:</h4>
	            <p class="footer__text text type--white">Universitat Pompeu Fabra</p>
	            <h4 class="footer__heading type--purple">Purpose:</h4>
	            <p class="footer__text text type--white">Use the CRISPR-A service provided by Translational Synthetic Biology research group</p>
	            <h4 class="footer__heading type--purple">Rights:</h4>
	            <p class="footer__text text type--white">you can access your data; request their rectification, deletion and portability; you may object to their processing and apply for their limitation</p>

	        </div>
	        <div class="footer__column">
	            <h4 class="footer__heading type--yellow">Contact:</h4>
	            <p class="footer__text text type--white">Marta Sanvicente Garcia:</p>
	            <a href="mailto:marta.sanvicente@upf.edu" class="footer__link footer__text text--bold type--white">marta.sanvicente at upf.edu</a>
	        </div>
	        <div class="footer__column footer__column--last">
	            <a href="documentation.html" class="footer__button button button--white">Need some help?</a>
	        </div>
	    </div>

	    <div class="footer-copyright footer__row">
	        <div class="footer__column">
	            <ul class="footer__list">
	                <li class="footer__item">
	                    <a href="../../TermsOfService.pdf" class="footer-copyright__text footer__link type--dark-grey">Legal Notice</a>
	                    <a href="../../PrivacyPolicy.pdf" class="footer-copyright__text footer__link footer__link--last type--dark-grey">Privacy Policy</a>
	                </li>
	            </ul>
	        </div>

	        <div class="footer__column">
	            <ul class="footer__list">1600
	                <li class="footer__item">
	                    <span class="footer-copyright__text footer-copyright__link type--dark-grey">CRISPR-A 2022</span>
	                </li>
	            </ul>
	        </div>

	        <div class="footer__column footer__column--last">
	            <ul class="footer__list">
	                <li class="footer__item">
	                    <span class="footer__text text type--white">Powered by <span class="footer__text text text--bold type--orange">Synbio</span></span>
	                </li>
	            </ul>
	        </div>
	    </div>

	</footer>

	
</body>

</html>
