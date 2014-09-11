<?php


function output($dir,$name){
        $context='<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
        <head>
        <title>TransPS</title>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
        <link href="style.css" rel="stylesheet" type="text/css" />
        <script type="text/javascript" src="form.js">
        
        </script>
        </head>
        <body>
        <div id="main">
        <div id="header">
        <div id="title">
        <!-- Name of tool -->
        <h1><a href="index.html"><strong>TransPS</strong>: Transcriptome Post Scaffolding. </a></h1>
        
        <!-- Authors -->
        <p>Mingming Liu, Zachary N Adelman, Liqing Zhang</p>
        
        <!-- Contact email -->
        <address>mingml@vt.edu</address>
        </div>
        </div>
        <div class="content">
        <div class="left">
        <div class="download">
        <a href="TransPS1.1.0.tar.gz"><img src="images/download_logo.png" alt="download" width="250"/></a>
        </div>
        <h2>Other data</h2>
        <ul>
        <li><a href="De_frontalis.tar.gz">De. frontalis</a></li>
        <li><a href="De_ponderosae.tar.gz">De. ponderosae</a></li>
        <li><a href="Di_citri.tar.gz">Di. citri</a></li>
        <li><a href="Ix_ricinus.tar.gz">Ix. ricinus</a></li>
        </ul>
        
        </div>
        <div class="right">
        <h2>Results</h2>
        <table id="table">
        <tr><td><p>TransPS generated the following output files,</p></td></tr>
        <tr><td>Accepted Contigs:</td><td>'."<a href=$dir/$name.accept>$name.accept </a>".'</td></tr>
        <tr><td>Scaffolded Contigs:</td><td>'."<a href=$dir/$name.scaffold>$name.scaffold </a>".'</td></tr>
        <tr><td>Unused Contigs:</td><td>'."<a href=$dir/$name.unused>$name.unused</a>".'</td></tr>
        <tr><td>Mapping:</td><td>'."<a href=$dir/$name.map>$name.map</a>".'</td></tr>
        </table>
        </div>
	<div class="footer">
	The website is under construction...
	</div>
        
        
        </div>
        
        </div>
        </body>
        </html>';
        
        /*$context1="<html><head><title>Results files</title></head>
         <body>
         <h2>Results files</h2>
         <ul>
         <li><a href=$dir/$name.map>$name:mapping </a></li>
         <li><a href=$dir/$name.accept>$name:accepted contigs </a></li>
         <li><a href=$dir/$name.scaffold>$name:scaffolding contigs</a></li>
         <li><a href=$dir/$name.unused>$name:unused contigs</a></li>
         <ul>
         </body>
         </html>";
         
         */
        
        echo $context;
    }
function generateRandomString($length = 6) {
        $characters = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
        $randomString = '';
        for ($i = 0; $i < $length; $i++) {
            $randomString .= $characters[rand(0, strlen($characters) - 1)];
        }
        return $randomString;
}

function myexec($cmd){
	exec($cmd,$output,$return);
//	var_dump($output); 
		if(count($output)<1){
			echo "Invalid arguments";
			return 1;
		}
		return 0;
//	echo $cmd;
}

function mycopy($orig,$dest){
//	echo "$orig<br>";
//	echo $dest;
	copy($orig,$dest) or die ("Unable to copy map file to $dest.");
}

$ttext = $_POST["seqs"];
$evalue = $_POST["evalue"];
$per = $_POST["per"];
$rate = $_POST["rate"];
$dist = $_POST["dist"];
$pwd = dirname(__FILE__);
$outdir="$pwd/results";


//echo $_FILES["seq_file"]["size"];
if(!empty($_FILES["seq_file"]["tmp_name"]) && !empty($_FILES["blastx_file"]["tmp_name"])){
	if($_FILES["seq_file"]["size"] < 20000000 &&  $_FILES["blastx_file"]["size"] < 20000000){
		if ($_FILES["seq_file"]["error"] > 0)
		{
			echo "Error: " . $_FILES["seq_file"]["error"] . "<br>";
			return;
		}
		if ($_FILES["blastx_file"]["error"] > 0)
		{
			echo "Error: " . $_FILES["blastx_file"]["error"] . "<br>";
			return;
		}
		$cmd = "perl TransPS1.1.0/transps.pl -t ".$_FILES["seq_file"]["tmp_name"]." -b ".$_FILES["blastx_file"]["tmp_name"]." --per $per --rate $rate --dist $dist --evalue $evalue";

//		echo $cmd;	
		if(myexec($cmd)==1){return;}

//		$fname = substr($_FILES["seq_file"]["name"],0,6);
		$name = $_FILES["seq_file"]["name"];
//		echo $dir;
        $str=generateRandomString();
        $dir = "$outdir/$str";
	if(!is_dir($dir)){
			$old=umask(0);
			mkdir("$dir",0757);
			umask($old);
		}
		mycopy($_FILES["seq_file"]["tmp_name"].".map","$dir/$name.map");
		mycopy($_FILES["seq_file"]["tmp_name"].".scaffold","$dir/$name.scaffold");
		mycopy($_FILES["seq_file"]["tmp_name"].".accept","$dir/$name.accept");
		mycopy($_FILES["seq_file"]["tmp_name"].".unused","$dir/$name.unused");
		output("https://bioinformatics.cs.vt.edu/zhanglab/transps/results/$str",$name);
	
		
	}
	
	else{
		echo "Invalide file!";
		return;
	}
}

else if(empty($_FILES["seq_file"]["tmp_name"]) && $ttext !== "" && !empty($_FILES["blastx_file"]["tmp_name"])){
	//convert text area input into fasta file and use uploaded blastx file
	if($_FILES["blastx_file"]["tmp_name"]<20000000){
		$name = $_FILES["blastx_file"]["name"].".fasta";
		$indir="/tmp";
		$fh = fopen("$indir/$name","w+");
		if($fh===false || fwrite($fh,$ttext)===false){
			echo "Unable to write file $indir/$name";
			return;
		}

		$cmd = "perl TransPS1.1.0/transps.pl -t $indir/$name -b ".$_FILES["blastx_file"]["tmp_name"]." --per $per --rate $rate --dist $dist --evalue $evalue";
		if(myexec($cmd)==1){return;}
		$str=generateRandomString();
		$dir="$outdir/$str";
		if(!is_dir($dir)){
			$old = umask(0);
			mkdir($dir,0777);
			umask($old);
		}
		mycopy("$indir/$name.map","$dir/$name.map");
		mycopy("$indir/$name.accept","$dir/$name.accept");
		mycopy("$indir/$name.scaffold","$dir/$name.scaffold");
		mycopy("$indir/$name.unused","$dir/$name.unused");
		output("https://bioinformatics.cs.vt.edu/zhanglab/transps/results/$str",$name);
		
	}
	else {
		echo "Invalid file";
	}
}

else if(empty($_FILES["seq_file"]["tmp_name"]) && empty($_FILES["blastx_file"]["tmp_name"]) && $ttext !== "" && isset($_POST['check'])){
	//convert text area input into fasta file and use example blastx file
	$name="example.fa";
	$indir="$pwd/results";
	$cmd = "perl TransPS1.1.0/transps.pl -t $indir/$name -b $indir/example.blastx --per $per --rate $rate --dist $dist --evalue $evalue";
//	echo $cmd;
	if(myexec($cmd)==1){return;}
	output("https://bioinformatics.cs.vt.edu/zhanglab/transps/results",$name);
}


else{
	echo "Invalide input!";
}


?>
