/*
The following provide the functionality for testing the VIP after updating using MS data
*/
var testVIPimg, testVIPlist;

function testVIP(eID){
  //$('#'+eID).animate({
  //  scrollTop: 0
  //},'slow');
  //$('html,body').scrollTop(0);
  //window.scrollTo(0, 0);

  // if the test has started already
  if($("#"+eID).html().includes("TESTING VIP")){
    return;
  }
  $("#"+eID).html("");
  testVIPhtmladd(eID,"<h3>TESTING VIP</h3>This might take some time, please be patient.<br>");

  testVIPgetData('static/testVIP/info.txt',testVIPload,eID);
  testVIPgetData('static/testVIP/img.txt',testVIPimgGet,eID);
  //All plots testing sequentially in testVIPimgGet function, since the async data retriving
}

function testVIPhtmladd(eID,message){
  $("#"+eID).html($("#"+eID).html()+message);
}

function testVIPgetData(url,callback){
  var args = Array.prototype.slice.call(arguments, 2);
  var xhr = new XMLHttpRequest();
  xhr.onload = function(){
    if(xhr.readyState === 4 && xhr.status === 200){
      callback.apply(xhr,args);
    }
  }
  xhr.open('GET', url, true);
  xhr.send(null);
}
function testVIPload(eID){
  var content = JSON.parse(decodeURIComponent(this.responseText));
  loadexe(content);
  testVIPhtmladd(eID,"Loading setting complete<br>");
}
function testVIPimgGet(eID){
  testVIPimg = JSON.parse(decodeURIComponent(this.responseText));
  testVIPhtmladd(eID,"Loading information complete<br>");
  testVIPlist = Object.keys(testVIPimg);
  testVIPnext(eID);
}

// The following testing all of them at once cause server 'Segmentation fault'
// Segmentation fault      (core dumped) cellxgene launch ...
function testVIPall(eID){
  for(const one in testVIPimg){
    var fn = window[one.replace('img','plot')];
    if(one.includes("DEG")){
      fn = window["DEGfind"];
    }
    if(typeof fn === "function") fn();
    testVIPhtmladd(eID,"Testing"+one.replace('img','')+" ...<br>");
    testVIPone(one,eID,0);
  }
}
function testVIPnext(eID){
  if(testVIPlist.length>0){
    var one = testVIPlist.pop();
    $("#"+one).remove();
    testVIPhtmladd(eID,"Testing "+one.replace('img','')+" ...<br>");
    var fn = window[one.replace('img','plot')];
    if(one.includes("DEG")){
      fn = window["DEGfind"];
    }
    if(typeof fn === "function") fn();
    testVIPone(one,eID,0);
  }else{
    testVIPhtmladd(eID,"=========== Testing complete! ===========<br>");
    testVIPhtmladd(eID,"<button onclick='$(\"#"+eID+"\").html(\"\");testVIP(\""+eID+"\")'>Retest</button>");
  }
}
function testVIPone(imgID,eID,nTimes){
  if(nTimes>60){//1 minutes
    testVIPhtmladd(eID,"<p style='color:red;'>Testing on "+imgID+" timeout</p>");
    testVIPnext(eID);
    return;
  }
  var info = $('#'+imgID).attr('src');
  if(info===undefined){
    setTimeout(testVIPone,1000,imgID,eID,nTimes+1);
  }else{
    if(info===testVIPimg[imgID]){
      testVIPhtmladd(eID,"Testing on "+imgID+" complete successfully<br>");
    }else if(info.length>100){
      testVIPhtmladd(eID,"<p style='color:salmon;'>Testing on "+imgID+" image changed</p>");
    }
    else{
      testVIPhtmladd(eID,"<p style='color:darkred;'>Testing on "+imgID+" failed</p>");
    }
    testVIPnext(eID);
  }
}
