// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "GSM437093", "1", "04/02/04 11:18:12" ], [ "2", "GSM437094", "2", "04/02/04 13:30:30" ], [ "3", "GSM437095", "3", "04/02/04 16:48:55" ], [ "4", "GSM437096", "4", "04/09/04 11:06:31" ], [ "5", "GSM437097", "5", "04/09/04 12:13:02" ], [ "6", "GSM437098", "6", "04/09/04 12:26:04" ], [ "7", "GSM437099", "7", "04/09/04 14:55:54" ], [ "8", "GSM437100", "8", "04/09/04 15:47:23" ], [ "9", "GSM437101", "9", "04/09/04 17:04:10" ], [ "10", "GSM437102", "10", "04/28/04 14:02:39" ], [ "11", "GSM437103", "11", "04/13/04 13:36:38" ], [ "12", "GSM437104", "12", "04/27/04 14:17:30" ], [ "13", "GSM437105", "13", "03/18/04 16:30:29" ], [ "14", "GSM437106", "14", "03/23/04 10:54:52" ], [ "15", "GSM437107", "15", "03/23/04 12:18:48" ], [ "16", "GSM437108", "16", "03/23/04 14:54:08" ], [ "17", "GSM437109", "17", "03/25/04 13:07:14" ], [ "18", "GSM437110", "18", "03/25/04 14:12:54" ], [ "19", "GSM437111", "19", "03/25/04 15:52:18" ], [ "20", "GSM437112", "20", "07/02/04 12:38:53" ], [ "21", "GSM437113", "21", "06/22/04 14:23:04" ], [ "22", "GSM437114", "22", "04/27/04 15:20:41" ], [ "23", "GSM437115", "23", "05/06/04 13:06:52" ], [ "24", "GSM437116", "24", "07/28/04 15:38:07" ], [ "25", "GSM437117", "25", "04/02/04 10:39:53" ], [ "26", "GSM437118", "26", "04/02/04 11:31:31" ], [ "27", "GSM437119", "27", "04/02/04 14:14:14" ], [ "28", "GSM437120", "28", "07/07/04 11:14:34" ], [ "29", "GSM437121", "29", "04/02/04 15:03:53" ], [ "30", "GSM437122", "30", "04/09/04 11:33:28" ], [ "31", "GSM437123", "31", "06/16/04 15:46:02" ], [ "32", "GSM437124", "32", "04/09/04 14:18:52" ], [ "33", "GSM437125", "33", "04/09/04 13:46:13" ], [ "34", "GSM437126", "34", "04/09/04 16:43:48" ], [ "35", "GSM437127", "35", "04/09/04 17:47:31" ], [ "36", "GSM437128", "36", "04/09/04 18:00:22" ], [ "37", "GSM437129", "37", "07/02/04 11:43:04" ], [ "38", "GSM437130", "38", "04/13/04 11:10:13" ], [ "39", "GSM437131", "39", "04/13/04 12:07:17" ], [ "40", "GSM437132", "40", "04/13/04 12:55:18" ], [ "41", "GSM437133", "41", "04/27/04 13:47:32" ], [ "42", "GSM437134", "42", "05/06/04 11:15:43" ], [ "43", "GSM437135", "43", "06/16/04 16:27:22" ], [ "44", "GSM437136", "44", "04/13/04 15:48:10" ], [ "45", "GSM437137", "45", "04/13/04 15:09:01" ], [ "46", "GSM437138", "46", "03/18/04 17:50:42" ], [ "47", "GSM437139", "47", "03/18/04 14:50:07" ], [ "48", "GSM437140", "48", "03/18/04 18:55:49" ], [ "49", "GSM437141", "49", "03/18/04 15:30:38" ], [ "50", "GSM437142", "50", "03/18/04 17:00:54" ], [ "51", "GSM437143", "51", "04/27/04 15:08:18" ], [ "52", "GSM437144", "52", "03/18/04 19:09:18" ], [ "53", "GSM437145", "53", "05/06/04 16:41:27" ], [ "54", "GSM437146", "54", "04/21/04 12:27:47" ], [ "55", "GSM437147", "55", "03/23/04 13:02:26" ], [ "56", "GSM437148", "56", "03/23/04 15:24:19" ], [ "57", "GSM437149", "57", "04/21/04 14:21:51" ], [ "58", "GSM437150", "58", "03/25/04 12:53:08" ], [ "59", "GSM437151", "59", "03/25/04 13:41:50" ], [ "60", "GSM437152", "60", "03/25/04 14:25:52" ], [ "61", "GSM437153", "61", "03/25/04 14:38:46" ], [ "62", "GSM437154", "62", "03/25/04 15:02:01" ], [ "63", "GSM437155", "63", "03/25/04 15:17:51" ], [ "64", "GSM437156", "64", "03/25/04 16:53:14" ], [ "65", "GSM437157", "65", "06/22/04 12:41:09" ], [ "66", "GSM437158", "66", "04/06/04 15:45:41" ], [ "67", "GSM437159", "67", "06/22/04 16:20:55" ], [ "68", "GSM437160", "68", "04/21/04 15:23:21" ], [ "69", "GSM437161", "69", "04/06/04 16:59:21" ], [ "70", "GSM437162", "70", "04/06/04 17:58:11" ], [ "71", "GSM437163", "71", "05/06/04 12:29:18" ], [ "72", "GSM437164", "72", "06/22/04 13:37:43" ], [ "73", "GSM437165", "73", "06/22/04 14:05:54" ], [ "74", "GSM437166", "74", "07/28/04 13:33:27" ], [ "75", "GSM437167", "75", "04/27/04 11:15:27" ], [ "76", "GSM437168", "76", "04/27/04 12:40:42" ], [ "77", "GSM437169", "77", "04/27/04 13:06:22" ], [ "78", "GSM437170", "78", "04/27/04 16:02:55" ], [ "79", "GSM437171", "79", "04/27/04 16:33:09" ], [ "80", "GSM437172", "80", "04/27/04 14:43:13" ], [ "81", "GSM437173", "81", "04/27/04 15:37:59" ], [ "82", "GSM437174", "82", "04/21/04 11:47:47" ], [ "83", "GSM437175", "83", "04/02/04 12:35:13" ], [ "84", "GSM437176", "84", "04/02/04 13:44:41" ], [ "85", "GSM437177", "85", "04/02/04 13:58:05" ], [ "86", "GSM437178", "86", "04/02/04 15:53:30" ], [ "87", "GSM437179", "87", "04/02/04 17:21:22" ], [ "88", "GSM437180", "88", "04/02/04 17:08:40" ], [ "89", "GSM437181", "89", "04/27/04 15:50:26" ], [ "90", "GSM437182", "90", "04/09/04 11:47:46" ], [ "91", "GSM437183", "91", "04/09/04 16:00:06" ], [ "92", "GSM437184", "92", "04/09/04 13:15:28" ], [ "93", "GSM437185", "93", "04/09/04 16:14:37" ], [ "94", "GSM437186", "94", "04/09/04 16:29:39" ], [ "95", "GSM437187", "95", "05/06/04 14:16:48" ], [ "96", "GSM437188", "96", "04/27/04 13:33:58" ], [ "97", "GSM437189", "97", "04/13/04 11:42:02" ], [ "98", "GSM437190", "98", "04/13/04 11:54:48" ], [ "99", "GSM437191", "99", "04/13/04 14:35:28" ], [ "100", "GSM437192", "100", "04/13/04 16:29:43" ], [ "101", "GSM437193", "101", "03/18/04 12:37:40" ], [ "102", "GSM437194", "102", "03/18/04 14:21:16" ], [ "103", "GSM437195", "103", "03/18/04 14:36:31" ], [ "104", "GSM437196", "104", "03/18/04 15:45:33" ], [ "105", "GSM437197", "105", "03/18/04 18:03:50" ], [ "106", "GSM437198", "106", "03/18/04 15:58:35" ], [ "107", "GSM437199", "107", "07/07/04 12:52:55" ], [ "108", "GSM437200", "108", "03/18/04 19:36:04" ], [ "109", "GSM437201", "109", "03/18/04 19:23:00" ], [ "110", "GSM437202", "110", "03/23/04 11:49:13" ], [ "111", "GSM437203", "111", "03/23/04 11:21:05" ], [ "112", "GSM437204", "112", "07/07/04 12:12:31" ], [ "113", "GSM437205", "113", "05/06/04 11:44:12" ], [ "114", "GSM437206", "114", "04/21/04 12:13:43" ], [ "115", "GSM437207", "115", "03/23/04 12:04:18" ], [ "116", "GSM437208", "116", "03/23/04 13:59:21" ], [ "117", "GSM437209", "117", "03/25/04 13:24:35" ], [ "118", "GSM437210", "118", "03/23/04 15:08:23" ], [ "119", "GSM437211", "119", "03/25/04 16:05:25" ], [ "120", "GSM437212", "120", "03/25/04 16:18:25" ], [ "121", "GSM437213", "121", "03/25/04 16:31:35" ], [ "122", "GSM437214", "122", "03/25/04 17:09:22" ], [ "123", "GSM437215", "123", "05/06/04 11:31:04" ], [ "124", "GSM437216", "124", "07/07/04 13:56:08" ], [ "125", "GSM437217", "125", "06/22/04 11:47:51" ], [ "126", "GSM437218", "126", "04/06/04 13:44:35" ], [ "127", "GSM437219", "127", "06/22/04 15:59:01" ], [ "128", "GSM437220", "128", "04/21/04 14:57:23" ], [ "129", "GSM437221", "129", "04/06/04 15:31:26" ], [ "130", "GSM437222", "130", "04/06/04 16:46:18" ], [ "131", "GSM437223", "131", "04/06/04 18:42:59" ], [ "132", "GSM437224", "132", "04/13/04 15:35:33" ], [ "133", "GSM437225", "133", "06/22/04 14:35:58" ], [ "134", "GSM437226", "134", "06/22/04 14:49:38" ], [ "135", "GSM437227", "135", "04/27/04 14:55:39" ], [ "136", "GSM437228", "136", "04/27/04 14:01:30" ], [ "137", "GSM437229", "137", "07/28/04 14:44:47" ], [ "138", "GSM437230", "138", "07/28/04 15:24:49" ], [ "139", "GSM437231", "139", "04/02/04 12:14:27" ], [ "140", "GSM437232", "140", "05/06/04 10:57:19" ], [ "141", "GSM437233", "141", "07/07/04 12:39:29" ], [ "142", "GSM437234", "142", "04/09/04 11:20:33" ], [ "143", "GSM437235", "143", "07/02/04 11:56:20" ], [ "144", "GSM437236", "144", "04/09/04 17:21:48" ], [ "145", "GSM437237", "145", "04/09/04 17:34:30" ], [ "146", "GSM437238", "146", "04/13/04 13:21:43" ], [ "147", "GSM437239", "147", "04/13/04 13:52:15" ], [ "148", "GSM437240", "148", "04/13/04 14:06:30" ], [ "149", "GSM437241", "149", "04/13/04 16:17:06" ], [ "150", "GSM437242", "150", "03/18/04 12:53:15" ], [ "151", "GSM437243", "151", "03/18/04 13:21:53" ], [ "152", "GSM437244", "152", "03/18/04 13:36:27" ], [ "153", "GSM437245", "153", "03/18/04 17:21:48" ], [ "154", "GSM437246", "154", "03/23/04 10:29:25" ], [ "155", "GSM437247", "155", "03/23/04 11:36:42" ], [ "156", "GSM437248", "156", "03/23/04 10:42:06" ], [ "157", "GSM437249", "157", "04/21/04 11:02:57" ], [ "158", "GSM437250", "158", "04/21/04 10:48:38" ], [ "159", "GSM437251", "159", "07/07/04 14:23:13" ], [ "160", "GSM437252", "160", "03/23/04 13:19:02" ], [ "161", "GSM437253", "161", "03/25/04 12:24:45" ], [ "162", "GSM437254", "162", "03/25/04 14:00:10" ], [ "163", "GSM437255", "163", "04/21/04 11:16:08" ], [ "164", "GSM437256", "164", "05/06/04 12:15:44" ], [ "165", "GSM437257", "165", "06/22/04 12:54:24" ], [ "166", "GSM437258", "166", "06/22/04 13:09:48" ], [ "167", "GSM437259", "167", "04/21/04 15:10:42" ], [ "168", "GSM437260", "168", "04/21/04 14:06:44" ], [ "169", "GSM437261", "169", "06/22/04 15:44:32" ], [ "170", "GSM437262", "170", "07/28/04 17:42:09" ], [ "171", "GSM437263", "171", "04/21/04 14:35:21" ], [ "172", "GSM437264", "172", "04/27/04 12:02:12" ], [ "173", "GSM437265", "173", "04/27/04 12:28:00" ], [ "174", "GSM437266", "174", "07/28/04 14:01:10" ], [ "175", "GSM437267", "175", "07/28/04 14:17:19" ], [ "176", "GSM437268", "176", "07/28/04 14:30:51" ], [ "177", "GSM437269", "177", "07/28/04 15:11:34" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
