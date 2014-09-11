<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" dir="ltr" lang="en-US" xml:lang="en">
    <head>
        <!--
         Base template (without user's data) checked by http://validator.w3.org : "This page is valid XHTML 1.0 Transitional"
         -->
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
        <meta http-equiv="X-UA-Compatible" content="IE=EmulateIE7" />
        <title>Help</title>
        
               
    </head>
    <body>
        
        <table width="405px" cellpadding="0" cellspacing="0" align="center">
            <tr><th bgcolor="#660000"><font size="4" color="#ffffff"> TransPS Online Help</font></th></tr>
            
            
            <tr><td bgcolor="#ffffff" align="center"><font size="3"></font></td></tr>
            <tr><td bgcolor="#ffffff">&nbsp;</td></tr>
            <tr><td bgcolor="#ffffff" style="text-align:justify;"><font size="2"></font>             
             <?php
                        $query = $_GET['page'];
                        if($query == 1){
				$info = "evalue: Expected threshold from BLASTX. If the evalue of a match is greater than a
					prespecified value, the contig will be droped into unused group, otherwise
					keep as accepted or used for scaffolding.";
                        }
                        else if($query == 2){
				$info = "overlapping percentage: If two contigs match with the same reference protein sequence and
					they overlapped greater than a given threshold, one of them will be moved to the unused group.";
                        }
                        else if($query == 3){
                        $info = "extention rate: In scaffolding case 2 (please see the paper), it is used for extending the scaffolding region.";
                        }
                        else if($query == 4){
                        $info = "scaffolding distance: During scaffolding, if the distance between two contigs is within the scaffolding distance, then join them together.";
                        }
            
                        
                        echo $info;
              ?>
            </td></tr>
        </table>
       
    </body>
</html>
