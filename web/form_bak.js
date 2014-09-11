 if (document.images) {
  arrow = new Image(); arrow.src = "images/arrow.png";
  arrow90 = new Image(); arrow90.src = "images/arrow90.png";
  }
  
  function hideShow(element){
    if((document.getElementById(element).style.display == 'none') || (document.getElementById(element).style.display == '')) {
      document.getElementById(element).style.display = 'block';
      if (document.images) {
        document.getElementById("arrow").src=arrow.src;
      }
    } else if(document.getElementById(element).style.display == 'block') {
      document.getElementById(element).style.display = 'none';
      if (document.images) {
        document.getElementById("arrow").src=arrow90.src;
      }
    }
  }
  
  function fill_example(){
  document.seq_form.seqs.value=">gi|417490306|gb|GACJ01016787.1| \n\
CCGGGTTCACTTGAATGGTGAAACAGTTCGGTTTGTCGGTCGAAACAAGCGAATAAG\nACAAGGCAATATT\
CGGATCAACAATCTGGACTGGGGTTTCACAAGGTTGGTACATTT\nGGCCTGACAACGGTGGTGATTCTTGT\
ACTTCTTGAGAACTGTATCCTTGATCTCCACC\nTGCTGGCTGATCAGCATTACCCAAGTTATTATGTGCGT\n\
>gi|417490306|gb|GACJ01016787.1| \n\
CCGGGTTCACTTGAATGGTGAAACAGTTCGGTTTGTCGGTCGAAACAAGCGAATAAG\nACAAGGCAATATT\
CGGATCAACAATCTGGACTGGGGTTTCACAAGGTTGGTACATTT\nGGCCTGACAACGGTGGTGATTCTTGT\
ACTTCTTGAGAACTGTATCCTTGATCTCCACC\nTGCTGGCTGATCAGCATTACCCAAGTTATTATGTGCGT\n";
  document.getElementById("check").checked=true;
  }
  
  function on_help(page){
  if(page==1)
  {
  window.open("help.php?page=1","Online Help",'height=250,width=480');
  }else if(page==2){
  window.open("help.php?page=2","Online Help",'height=250,width=480');
  }else if(page==3){
  window.open("help.php?page=3","Online Help",'height=250,width=480');
  }else{
  window.open("help.php?page=4","Online Help",'height=250,width=480');
  }
  
  }
  
  function validateForm(){
  	  var ttext = document.forms["seq_form"]["seqs"].value;
  	  var tfile = document.forms["seq_form"]["seq_file"].value;
  	  var bfile = document.forms["seq_form"]["blastx_file"].value;
  	  var evalue = document.forms["seq_form"]["evalue"].value;
  	  var per = document.forms["seq_form"]["per"].value;
  	  var rate = document.forms["seq_form"]["rate"].value;
  	  var dist = document.forms["seq_form"]["dist"].value;
  	  if((ttext==null||ttext==""||/^\s+$/g.test(ttext)) && (tfile==null||tfile=="")){
  	  	  alert("Please input/upload transcriptome file!");
  	  	  return false;
  	  }
  	  if(evalue==""||per==""||rate==""||dist==""){
  	  	  alert("Please set the arguments!");
  	  	  return false;
  	  }
  	  if((bfile==null||bfile=="") && document.getElementById("check").checked==false){
  	  	  alert("Please upload BLASTX output file!");
  	  	  return false;
  	  }
  }
  
  
