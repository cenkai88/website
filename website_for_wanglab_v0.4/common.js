//检查输入的氨基酸比例是否在0到1之间
function checkaaratio (ratio) {
 	if (ratio > 1 || ratio < 0) {
 		return false
 	}
 }
 //搜索蛋白组，长度，信号肽，氨基酸比例
function search (proteome, len, signalp, aa, ratio, tmnumber) {
	var tmp = new Array()
	if (signalp == 'Y') {
		for (var i = 0; i < proteome.length; i++ ){
			var n = (proteome[i].seq.split(aa).length - 1 ) / proteome[i].seq.length
			if (proteome[i].seq.length <= len && proteome[i].signal == 'Y' && n >= ratio && tmjudge(proteome[i].tm, tmnumber)){
				tmp.push(proteome[i])
			}
		}
	}
	else {
		for (var i = 0; i < proteome.length; i++ ){
			var n = (proteome[i].seq.split(aa).length - 1 ) / proteome[i].seq.length
			if (proteome[i].seq.length <= len && n >= ratio && tmjudge(proteome[i].tm, tmnumber)){
				tmp.push(proteome[i])
			}
		}
	}
	return tmp
}
//输出结果
function outputpros (result){
	outputWindow.document.write("---------------------------------------------------------------------------------------" + "<br>")
	outputWindow.document.write('>' + result.name + "<br>")
	x = Math.floor((result.seq.length-1) / 60)
	for (var i = 0; i < x + 1; i++ ){
		outputWindow.document.write(result.seq.substring(i*60, i*60+60) + "<br>")
	}
}
//物种名映射
function species (name){
	if (name == 'BBA'){
		return BBA
	}
	if (name == 'MAA'){
		return MAA
	}
	if (name == 'CCM'){
		return CCM
	}
	if (name == 'MAC'){
		return MAC
	}
	else{
		alert('尚无此物种信息，或输入格式不对，例：MAA_00001')
	}
}
//跨膜域
function tmjudge (tm, number){
	if (number == 'none'){
		return true
	}
	if (number =='有'){
		if (tm > 0) {
			return true
		}
		else {
			return false
		}
	}
	if (tm == number){
		return true
	}
}

//查询蛋白序列
function showseqinfo (spnumber){
	sp = spnumber.split("_")[0]
	number = spnumber.split("_")[1]
	proteome = species(sp)
	for (var i = 0; i < proteome.length ; i++){
		if (proteome[i].name.substr(-5) == number) {
			return outputproinfo(proteome[i])
		}
	}
}

//输出蛋白序列信息
function outputproinfo (result){
	tmp = "---------------------------------------------------------------------------------------" + "<br>"
	tmp = tmp + "基因名：" + result.name + "<br>"
	tmp = tmp + "长度：" + result.length + "<br>"
	x = Math.floor((result.seq.length-1) / 60)
	tmp = tmp + '序列：'
	for (var i = 0; i < x + 1; i++ ){
		tmp = tmp + result.seq.substring(i*60, i*60+60) + "<br>"
	}
	tmp = tmp + '是否有信号肽(signal peptide): ' + result.signal + "<br>"
	tmp = tmp + '跨膜域(transmembrane helices)数量: ' + result.tm + "<br>"
	tmp = tmp + '预测分子量(molecular weight): ' + mw(result) + '<br>'
	return tmp
}

//查询dna序列
function showseqinfodna (spnumber){
	sp = spnumber.split("_")[0]
	number = spnumber.split("_")[1]
	proteome = species(sp)
	for (var i = 0; i < proteome.length ; i++){
		if (proteome[i].name.substr(-5) == number) {
			return outputdnainfo(proteome[i])
		}
	}
}

//输出dna序列信息
function outputdnainfo (result){
	tmp = "----------------------------------------------------------------------------------------------------" + "<br>"
	tmp = tmp + "基因名：" + result.name + "<br>"
	tmp = tmp + '长度：' + result.dna.length + "<br>"
	tmp = tmp + 'cDNA长度：' + (parseInt(result.dna.length) - parseInt(intronlen(result))) + "<br>"
	x = Math.floor((result.dna.length-1) / 60)
	y = result.intron.split(',')
	tmp = tmp + '序列：'
	tmp2 = result.dna
	if (y.length > 1) {
		for (k = 0; k < y.length; k++){
			if (k%2 == 0){
				tmp2 = tmp2.substr(0, parseInt(y[k]) + parseInt(23*k/2)  ) + "<font color=red>" + tmp2.substr(parseInt(y[k])+ parseInt(23*k/2) )
			}
			else{
				tmp2 = tmp2.substr(0,parseInt(y[k]) + parseInt(4.5) + parseInt(11.5*k) ) + "</font>" + tmp2.substr(parseInt(y[k]) + parseInt(4.5) + parseInt(11.5*k))
			}
		}
	}
	if (tmp2.substr(0,3)=="ATG"){
		tmp2 = "<font color=green>" + tmp2.substr(0,3) + "</font>" + tmp2.substr(3)
		tmp2 =  tmp2.substr(0,tmp2.length-3) + "<font color=purple>"  + tmp2.substr(tmp2.length-3) +"</font>"
	}
	else {
		tmp2 = "<font color=purple>" + tmp2.substr(0,3) + "</font>" + tmp2.substr(3)
		tmp2 =  tmp2.substr(0,tmp2.length-3) + "<font color=green>"  + tmp2.substr(tmp2.length-3) +"</font>"
	}
	tmp2 = tmp2 +'<br>'
	s=1
	j=0
	while (tmp2[j]!='b'){
		if (s%60 ==0){
			tmp2 = tmp2.substr(0,j) + "<br>" + tmp2.substr(j)
			j = j + 4
		}
		if (tmp2[j]=='A' || tmp2[j]=='T' || tmp2[j]=='C' || tmp2[j]=='G'){
			s = s + 1
		}
		j = j + 1
	}
	tmp = tmp + tmp2  + '<br>'
	tmp = tmp + 'Contig：' + result.contig + "<br>"
	tmp = tmp + 'Contig前一个基因距离：' + result.former + "<br>"
	tmp = tmp + 'Contig后一个基因距离：' + result.next + "<br>"
	tmp = tmp + '内含子数目：' + Math.floor(result.intron.split(',').length / 2) + "<br>"
	tmp = tmp + "<br>如果Contig为none的话，说明本基因位于一个注释质量不高的Contig上，无法上传到NCBI" 
	return tmp
}

//统计氨基酸组成
function aa_ratio (protein){
	a = protein.seq.split('A').length - 1 
	v = protein.seq.split('V').length - 1 
	l = protein.seq.split('L').length - 1 
	i = protein.seq.split('I').length - 1 
	f = protein.seq.split('F').length - 1  
	w = protein.seq.split('W').length - 1 
	m = protein.seq.split('M').length - 1 
	p = protein.seq.split('P').length - 1 
	g = protein.seq.split('G').length - 1 
	s = protein.seq.split('S').length - 1 
	t = protein.seq.split('T').length - 1 
	c = protein.seq.split('C').length - 1 
	y = protein.seq.split('Y').length - 1 
	n = protein.seq.split('N').length - 1 
	q = protein.seq.split('Q').length - 1 
	h = protein.seq.split('H').length - 1 
	k = protein.seq.split('K').length - 1 
	r = protein.seq.split('R').length - 1 
	d = protein.seq.split('D').length - 1 
	e = protein.seq.split('E').length - 1 
	tmp = [a, v, l ,i, f, w, m, p, g, s, t, c, y, n, q, h, k, r, d, e]
	return tmp
}

//计算分子量
function mw(protein){
	weight = [71.0788, 99.1326, 113.1595, 113.1595, 147.1766, 186.2133, 131.1986, 97.1167, 57.0520, 87.0782, 101.1051, 103.1448, 163.1760, 114.1039, 128.1308, 137.1412, 128.1742, 156.1876, 115.0886, 129.1155]
	sum = 0
	for (var i = 0; i <20; i++){
		sum = sum + aa_ratio(protein)[i] * weight[i]
	}
	return sum/1000
}

// 计算内含子长度
function intronlen(protein){
	x = protein.intron.split(',')
	s = 0
	if (x.length == 1){
		return 0
	}
	for (var i=0; i<x.length/2; i++){
		s = parseInt(s) + parseInt(x[i*2+1]) - parseInt(x[i*2]) - parseInt(1)
	}
	return s
}



//Copyright (C) 2000 Paul Stothard
function addCodon (residueArray,iValue)	{
	codonString = codon[iValue];
	codonProb = fraction[iValue];
	for (var k = 0; k < codonString.length; k++)	{
		if (codonString.charAt(k) == "g")	{
			residueArray[k*4] = residueArray[k*4] + codonProb;
		}
		else	{
			if (codonString.charAt(k) == "a")	{
				residueArray[k*4+1] = residueArray[k*4+1] + codonProb;
			}
			else	{
				if (codonString.charAt(k) == "t")	{
					residueArray[k*4+2] = residueArray[k*4+2] + codonProb;
				}
				else	{
					if (codonString.charAt(k) == "c")	{
						residueArray[k*4+3] = residueArray[k*4+3] + codonProb;
					}
				}
			}
		}
	}
	return residueArray;
}
//在表达序列中间添加X，以便后续翻译
function addX (expressionString)	{	//adds X's to regular expression so that split function can be used. Expression needs both slashes.
	var expChar = "";
	var j = 1;
	while (expChar != "/")	{
		var tempMatchExp = expressionString.substring(j,j+2);               
		if (tempMatchExp.search(/[a-z]{2}/i) != -1)	{
			expressionString = expressionString.substring(0,j+1) + "X" + expressionString.substring(j+1,expressionString.length);
			j = j + 2;
		}
		else	{
			j = j + 1;
		}
		if (expChar == "[")	{	//handles square brackets in regular expressions (don't want 'X' within them).
			expressionString = expressionString.substring(0,j - 1) + "X" + expressionString.substring(j - 1,expressionString.length);	//places 'X' before square bracket
			while (expChar != "]")	{	
				j = j + 1;
				expChar = expressionString.charAt(j);
			}
			expressionString = expressionString.substring(0,j + 1) + "X" + expressionString.substring(j + 1,expressionString.length);	//places 'X' after square bracket.
			j = j + 1;	
		}
		expChar = expressionString.charAt(j);
	}
	expressionString = expressionString.replace(/XX/g,"X");	//remove double X's when square brackets next to eachother.
	expressionString = expressionString.replace(/\/X/,"/"); //remove 'X' from start of expression.
	expressionString = expressionString.replace(/X\//,"/"); //remove 'X' from end of expression.
	return expressionString;
}
function amAcidProb () {
	rA = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rC = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rD = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rE = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rF = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rG = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rH = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rI = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rK = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rL = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rM = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rN = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rP = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rQ = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rR = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rS = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rT = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rV = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rW = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rY = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	rZ = new Array(0,0,0,0,0,0,0,0,0,0,0,0);
	var arrayOfAcceptedAminoAcids = ["Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu","Met","Asn","Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr","End"];
	var arrayOfAminoAcidProb = [rA,rC,rD,rE,rF,rG,rH,rI,rK,rL,rM,rN,rP,rQ,rR,rS,rT,rV,rW,rY,rZ];
	for (var i = 0; i < amAcid.length; i++)	{
		for (var j = 0; j < arrayOfAcceptedAminoAcids.length; j++)	{
			if (amAcid[i] == arrayOfAcceptedAminoAcids[j])	{
				arrayOfAminoAcidProb[j] = addCodon(arrayOfAminoAcidProb[j],i);
				break;
			}
		}
	}
	return true;
}
function checkAlign (arrayOfTitles,arrayOfSequences)	{
	var lengthOfAlign = arrayOfSequences[0].length;
	if (arrayOfSequences.length < 2)	{
		alert ("Please enter an alignment consisting of at least two sequences.");
		return false;
	}
	for (var i = 0; i < arrayOfTitles.length; i++)	{
		if ((arrayOfTitles[i].search(/\S/) == -1) || (arrayOfSequences[i].search(/\S/) == -1) || (arrayOfSequences[i].length != lengthOfAlign))	{
			alert ("There is a problem with the alignment format.");
			return false;
		}
	}
	return true;
}
function checkCodonTable (codonTable)	{
	if ((codonTable.search(/AmAcid/) == -1) || (codonTable.search(/Codon/) == -1) || (codonTable.search(/Number/) == -1) || (codonTable.search(/\/1000/) == -1) || (codonTable.search(/Fraction\s*\.\./) == -1)) {
		alert ("The codon table has been entered incorrectly.");
		return false;
	}
	return true;
}
//检查输入框是否为空
function checkFormElement (formElement)	{
	if (formElement.value.search(/\S/) == -1)	{
		alert ("The form has not been filled out completely.");
		return false;
	}
	return true;
}
//检查自定义密码子输入是否符合格式
function checkGeneticCode (arrayOfPatterns)	{
	var z = 0;
	var codon = "";
	var oneMatch = false;
	var testSequence = "gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaagtagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagtaatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcacctgctacttctcccgccacctccc";
	while (z < arrayOfPatterns.length)	{
		if (arrayOfPatterns[z].search(/^\s*\/[a-zA-Z\|\[\]]+\/=[a-zA-Z\*]/) == -1)	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		if (moreExpressionCheck(arrayOfPatterns[z]) == false)	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		z = z + 1;
	}
	geneticCodeMatchResult = new Array (arrayOfPatterns.length);	//can be used globally for translation-related functions.
	geneticCodeMatchExp = new Array (arrayOfPatterns.length);	//can be used globally for translation-related functions.
	for (var j = 0; j < arrayOfPatterns.length; j++)	{
			geneticCodeMatchExp[j] = eval(arrayOfPatterns[j].match(/\/.+\//) + "gi");
			geneticCodeMatchResult[j] = (arrayOfPatterns[j].match(/=[a-zA-Z\*]/)).toString();
			geneticCodeMatchResult[j] = geneticCodeMatchResult[j].replace(/=/g,"");
	}
	for (var i = 0; i <= testSequence.length - 3; i = i + 3)	{
		codon = testSequence.substring(i,(i+3));
		for (var j = 0; j < geneticCodeMatchExp.length; j++)	{
			if (codon.search(geneticCodeMatchExp[j]) != -1)	{
				if (oneMatch == true)	{
					alert('More than one amino acid is coded by a single codon.');
					return false;
				}
				oneMatch = true;
			}
		}
		if (oneMatch == false)	{
			alert('The genetic code expressions are missing a codon.');
			return false;
		}
		oneMatch = false;
	}
	return true;
}
function checkGroupInput (arrayOfPatterns)	{
	var z = 0;
	while (z < arrayOfPatterns.length)	{
		if (arrayOfPatterns[z].search(/[^acdefghiklmnpqrstvwyz]/i) != -1)	{	
			alert ("One or more groups have been entered incorrectly.");
			return false;
		}
		z = z + 1;
	}
	for (var i = 0; i < arrayOfPatterns.length; i++)	{
		var re = new RegExp ("[" + arrayOfPatterns[i] + "]","gi");
		for (var j = i + 1; j < arrayOfPatterns.length; j++)	{
			if (arrayOfPatterns[j].search(re) != -1)	{
				alert('The same amino acid is in more than one similarity group.');
				return false;
			}
		}
	}
	return true;
}
function checkIdent (columnSeq,aminoAcid)	{
	var re = new RegExp (aminoAcid,"gi");
	return (columnSeq.match(re)).length;
}
function checkPatternInput (arrayOfPatterns)	{
	var z = 0;
	while (z < arrayOfPatterns.length)	{
		if ((arrayOfPatterns[z].search(/^\s*\/[GATCNgatcn\[\]]+\/\s+\([^\/]+\)/) == -1) || (arrayOfPatterns[z].search(/\(\)/) != -1))	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		if (moreExpressionCheck(arrayOfPatterns[z]) == false)	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		z = z + 1;
	}
	return true;
}
function checkProtPatternInput (arrayOfPatterns)	{
	var z = 0;
	while (z < arrayOfPatterns.length)	{
		if ((arrayOfPatterns[z].search(/^\s*\/[ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\[\]]+\/\s+\([^\/]+\)/) == -1) || (arrayOfPatterns[z].search(/\(\)/) != -1))	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		if (moreExpressionCheck(arrayOfPatterns[z]) == false)	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		z = z + 1;
	}
	return true;
}
//检查自定义酶切位点的输入是否符合格式
function checkRestPatterns (arrayOfPatterns)	{
	var z = 0;
	while (z < arrayOfPatterns.length)	{
		if (arrayOfPatterns[z].search(/^\s*\/[GATCNgatcn\[\]]+\/\s+\([^\/]+\)\d+/) == -1)	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		if (moreExpressionCheck(arrayOfPatterns[z]) == false)	{
			alert ("One or more patterns have been entered incorrectly.");
			return false;
		}
		z = z + 1;
	}
	return true;
}
function checkSim (columnSeq,aminoAcid,arrayOfGroups)	{
	var result = 1;
	var re = new RegExp (aminoAcid,"gi");
	for (var m = 0; m < arrayOfGroups.length; m++)	{
		if (arrayOfGroups[m].search(re) != -1)	{
			var re = new RegExp ("[" + arrayOfGroups[m] + "]","gi");
			result = (columnSeq.match(re)).length;
			break;
		}
	}
	return result;
}
function closeForm ()	{
	outputWindow.document.write ("</FORM>");
	return true;
}
function closeTextArea ()	{
	outputWindow.document.write ('</TEXTAREA>');
	return true;
}
function closeWindow ()	{
	outputWindow.document.write ("</body></html>");
	outputWindow.status = 'Done.';
	outputWindow.document.close();
	return true;
}
function complement (dnaSequence)	{
	dnaSequence = dnaSequence.replace(/g/g,"x");
	dnaSequence = dnaSequence.replace(/G/g,"X");
	dnaSequence = dnaSequence.replace(/a/g,"y");
	dnaSequence = dnaSequence.replace(/A/g,"Y");
	dnaSequence = dnaSequence.replace(/t/g,"a");
	dnaSequence = dnaSequence.replace(/T/g,"A");
	dnaSequence = dnaSequence.replace(/c/g,"g");
	dnaSequence = dnaSequence.replace(/C/g,"G");
	dnaSequence = dnaSequence.replace(/x/g,"c");
	dnaSequence = dnaSequence.replace(/X/g,"C");
	dnaSequence = dnaSequence.replace(/y/g,"t");
	dnaSequence = dnaSequence.replace(/Y/g,"T");		
	return dnaSequence;
}
function complementExp (expressionString)	{	//complements the expression. Expression needs both slashes. Used in conjunction with the reverseExp function.
	var param = expressionString.match(/[a-z]+$/);
	var nonParam = expressionString.match(/\/.+\//);
	nonParam = nonParam.toString();
	nonParam = nonParam.replace(/g/gi,"1");
	nonParam = nonParam.replace(/a/gi,"2");
	nonParam = nonParam.replace(/t/gi,"a");
	nonParam = nonParam.replace(/c/gi,"g");
	nonParam = nonParam.replace(/1/gi,"c");
	nonParam = nonParam.replace(/2/gi,"t");
	nonParam = nonParam + param;
	return nonParam;
}
function earlyCheckAlign (alignArray)	{
	if (alignArray.length < 3)	{
		alert ("There is a problem with the alignment format.");
		return false;
	}
	for (var i = 1; i < alignArray.length; i++)	{
		if (alignArray[i].search(/[^\s]+\s/) == -1)	{
			alert ("There is a problem with the alignment format.");
			return false;
		}
	}
	return true;
}
function filterAlignSeq (alignSeq)	{
	alignSeq = alignSeq.replace(/[^acdefghiklmnpqrstvwyx\.\-]/gi,"");
	return alignSeq;
}
function filterDna (dnaSequence)	{
	return (dnaSequence.replace (/[^gatcn]/gi,"")).toLowerCase();
}
function filterDnaLight (dnaSequence)	{
	return (dnaSequence.replace (/[\s\d]/gi,"")).toLowerCase();
}
//序列转为小写
function filterDnaRegExp (expressionString)	{	
	expressionString = expressionString.toLowerCase();
	return expressionString;
}
//将DNA序列中gatcn外的字母查找并去掉
function filterDnaSaveCase (dnaSequence)	{
	return (dnaSequence.replace (/[^gatcn]/gi,""));
}
//去掉fasta标题前面的空白字符和>或<
function filterFastaTitle (sequenceTitle)	{
	sequenceTitle = sequenceTitle.replace(/\s{2,}/g," ");
	return sequenceTitle.replace (/[\<\>]/gi,"");
}
function filterProtein (proteinSequence)	{
	return (proteinSequence.replace (/[^acdefghiklmnpqrstvwyz\*]/gi,"")).toUpperCase();
}
function filterProteinLight (proteinSequence)	{
	return (proteinSequence.replace (/[^A-Za-z]/gi,"")).toUpperCase();
}
function filterProteinSaveCase (proteinSequence)	{
	return (proteinSequence.replace (/[^acdefghiklmnpqrstvwyz\*]/gi,""));
}
function formatTable (codonTable)	{
	amAcid = new Array;
	codon = new Array;
	number = new Array;
	perThou = new Array;
	fraction = new Array;
	tableArray = new Array;
	for (var i = 0; i < 64; i++)	{
		amAcid[i] = "1";
		codon[i] = "1";
		number[i] = "X";
		perThou[i] = "X";
		fraction[i] = "X";
	}
	var perThouSum = 0;
	var i = 0;
	var codonString = "";
	var acceptedAminoAcids = "Ala Cys Asp Glu Phe Gly His Ile Lys Leu Met Asn Pro Gln Arg Ser Thr Val Trp Tyr End";
	var arrayOfAcceptedAminoAcids = ["Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu","Met","Asn","Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr","End"];
	var arrayOfFractionValues = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
	codonCounter = 0;
	codonTableGood = true;
	codonTable = codonTable.replace(/[^\.]*\.\./,"");
	tableArray = codonTable.split(/[\f\n\r]/);
	for (var j = 0; j < tableArray.length; j++)	{
		if (tableArray[j].search(/[^\s]/) != -1)	{
			tempString = tableArray[j];
			if (tempString.search(/\S+/) != -1)	{
				amAcid[i] = (tempString.match(/\S+/)).toString();
				amAcid[i] = (amAcid[i].substring(0,1)).toUpperCase() + (amAcid[i].substring(1,3)).toLowerCase();
				tempString = tempString.replace(/\S+/,"");
			}
			if (tempString.search(/\S+/) != -1)	{
				codon[i] = tempString.match(/\S+/);
				codon[i] = (codon[i].toString()).toLowerCase();
				tempString = tempString.replace(/\S+/,"");
			}
			if (tempString.search(/\S+/) != -1)	{
				number[i] = parseFloat((tempString.match(/\S+/)).toString());
				tempString = tempString.replace(/\S+/,"");
			}
			if (tempString.search(/\S+/) != -1)	{
				perThou[i] = parseFloat((tempString.match(/\S+/)).toString());
				tempString = tempString.replace(/\S+/,"");
			}
			if (tempString.search(/\S+/) != -1)	{
				fraction[i] = parseFloat((tempString.match(/\S+/)).toString());
			}
			if ((amAcid[i].search(/[A-Z][a-z][a-z]/) == -1) || (codon[i].search(/[a-z][a-z][a-z]/) == -1) || (((number[i]).toString()).search(/\d/) == -1) || (((perThou[i]).toString()).search(/\d/) == -1) || (((fraction[i]).toString()).search(/\d/) == -1))	{
				codonTableGood = false;
				break;
			}
			var re = new RegExp (codon[i]);
			if (codonString.search(re) != -1)	{
				alert('The same codon appears twice in the codon table.');
				codonTableGood = false;
				break;
			}
			var re = new RegExp (amAcid[i]);
			if (acceptedAminoAcids.search(re) == -1)	{
				alert('There is an unrecognized amino acid in the codon table.');
				codonTableGood = false;
				break;
			}
			for (var k = 0; k < arrayOfAcceptedAminoAcids.length; k++)	{
				if (amAcid[i] == arrayOfAcceptedAminoAcids[k])	{
					arrayOfFractionValues[k] = arrayOfFractionValues[k] + fraction[i];
					break;
				}
			}
			perThouSum = perThouSum + perThou[i];
			codonString = codonString + codon[i] + " ";
			i = i + 1;
			codonCounter = codonCounter + 1;
		}
	}
	for (var k = 0; k < arrayOfFractionValues.length; k++)	{
		if (Math.abs(arrayOfFractionValues[k] - 1) > 0.03)	{
			codonTableGood = false;
			break;
		}
	}
	if (Math.abs(perThouSum - 1000) > 0.32)	{
		codonTableGood = false;
	}
	return true;
}
function formatTableForCode (codonTable)	{
	amAcid = new Array;
	codon = new Array;
	number = new Array;
	perThou = new Array;
	fraction = new Array;
	tableArray = new Array;
	for (var i = 0; i < 64; i++)	{
		amAcid[i] = "1";
		codon[i] = "1";
		number[i] = 0;
		perThou[i] = 0;
		fraction[i] = 0;
	}
	var i = 0;
	var codonString = "";
	var acceptedAminoAcids = "Ala Cys Asp Glu Phe Gly His Ile Lys Leu Met Asn Pro Gln Arg Ser Thr Val Trp Tyr End";
	codonCounter = 0;
	codonTableGood = true;
	codonTable = codonTable.replace(/[^\.]*\.\./,"");
	tableArray = codonTable.split(/[\f\n\r]/);
	for (var j = 0; j < tableArray.length; j++)	{
		if (tableArray[j].search(/[^\s]/) != -1)	{
			tempString = tableArray[j];
			if (tempString.search(/\S+/) != -1)	{
				amAcid[i] = (tempString.match(/\S+/)).toString();
				amAcid[i] = (amAcid[i].substring(0,1)).toUpperCase() + (amAcid[i].substring(1,3)).toLowerCase();
				tempString = tempString.replace(/\S+/,"");
			}
			if (tempString.search(/\S+/) != -1)	{
				codon[i] = tempString.match(/\S+/);
				codon[i] = (codon[i].toString()).toLowerCase();
				tempString = tempString.replace(/\S+/,"");
			}
			if ((amAcid[i].search(/[A-Z][a-z][a-z]/) == -1) || (codon[i].search(/[a-z][a-z][a-z]/) == -1))	{
				codonTableGood = false;
				break;
			}
			var re = new RegExp (codon[i]);
			if (codonString.search(re) != -1)	{
				alert('The same codon appears twice in the codon table.');
				codonTableGood = false;
				break;
			}
			var re = new RegExp (amAcid[i]);
			if (acceptedAminoAcids.search(re) == -1)	{
				alert('There is an unrecognized amino acid in the codon table.');
				codonTableGood = false;
				break;
			}
			codonString = codonString + codon[i] + " ";
			i = i + 1;
			codonCounter = codonCounter + 1;
		}
	}
	return true;
}
function groupNum(text,tabIn,groupSize,basePerLine,startBase,stopBase,sepChar)	{
	var i = startBase;
	var k = 0;
	var lineOfText = "";
	groupSize = parseFloat(groupSize);
	basePerLine = parseFloat(basePerLine);
	while (i < stopBase)	{
		lineOfText = rightNum(i + 1,lineOfText,8,tabIn);
		for (var j = 1; j <= (basePerLine/groupSize); j++)	{//makes a group each loop
			while (k < groupSize)	{
				lineOfText = lineOfText + text.charAt(k + i);
				k = k + 1;
			}
			i = i + groupSize;
			k = 0;
			lineOfText = lineOfText + sepChar;
		}
		outputWindow.document.write (lineOfText + "<br>\n");
		lineOfText = "";
	}
	return true;
}
function groupNumALine(text,tabIn,groupSize,startBase,sepChar)	{
	var i = 0;
	var k = 0;
	var lineOfText = "";
	var basePerLine = text.length;
	groupSize = parseFloat(groupSize);
	while (i < basePerLine)	{
		lineOfText = rightNum(startBase + 1,lineOfText,8,tabIn);
		for (var j = 1; j <= (basePerLine/groupSize); j++)	{//makes a group each loop
			while (k < groupSize)	{
				lineOfText = lineOfText + text.charAt(k + i);
				k = k + 1;
			}
			i = i + groupSize;
			k = 0;
			lineOfText = lineOfText + sepChar;
		}
		outputWindow.document.write (lineOfText + "<br>\n");
		lineOfText = "";
	}
	return true;
}
//用[gatcn]替换[gatc]
function handleN (expressionString)	{	//replaces [gatc] with [gatcn].  Must use toMinimum function first.
	expressionString = expressionString.replace(/\[gatc\]/g,"[gatcn]");
	return expressionString;
}
function makeTag (alignChar,amAcids,amColors,frontOrBack)	{
	var colorForTag = "";
	for (var k = 0; k < amAcids.length; k++)	{
		if (alignChar.toUpperCase() == amAcids[k].toUpperCase())	{
			colorForTag = amColors[k];
			break;
		}
	}
	var amColorsArray = amColors[k].split(/\,\s/);
	if (frontOrBack == 'front')	{
		return '<font color="' + amColorsArray[0] + '">' + alignChar + '</font>';
	}
	else	{
		return '<span style="background-color:' + amColorsArray[0] + '"><font color="' + amColorsArray[1] + '">' + alignChar + '</font></span>';
	}
}
//检查是否有[字母[，]字母],[],/字母],[字母/,||,/|,|/,[任意],>,<,如有，则返回F
function moreExpressionCheck (stringToCheck)	{
	if ((stringToCheck.search(/\[[A-Za-z\|]*\[/) != -1) || (stringToCheck.search(/\][A-Za-z\|]*\]/) != -1) || (stringToCheck.search(/\[\]/) != -1) || (stringToCheck.search(/\/[A-Za-z\|]*\]/) != -1) || (stringToCheck.search(/\[[A-Za-z\|]*\//) != -1) || (stringToCheck.search(/\|\|/) != -1) || (stringToCheck.search(/\/\|/) != -1) || (stringToCheck.search(/\|\//) != -1) || (stringToCheck.search(/\[.\]/) != -1) || (stringToCheck.search(/\</) != -1) || (stringToCheck.search(/\>/) != -1))	{
		return false;
	}
	return true;
}
function openForm ()	{
	outputWindow.document.write ('<FORM NAME="output">');
	return true;
}
function openTextArea ()	{
	outputWindow.document.write ('<br><TEXTAREA NAME="output" ROWS="6" COLS="60" WRAP=SOFT>');
	return true;
}
function openWindow (title) {
	outputWindow=window.open("","displayWindow","resizable=yes,menubar=yes,toolbar=yes,height=800,width=1000,scrollbars=yes,status=yes,offscreenBuffering=false");
	outputWindow.focus();
	outputWindow.document.write ("<HTML><HEAD><TITLE>Output Window</TITLE>");
	outputWindow.document.write ("</HEAD>");
	outputWindow.document.write ('<body text="#000000" bgcolor="#FFFFFF" link="#0000FF" vlink="#551A8B" alink="#0000FF">');
	outputWindow.document.write ('<font size=+1>' + title + '</font><br>');
	outputWindow.status = 'Please Wait.';
	return true;
}
function palinCheck (expressionString)	{	//checks whether expression is palindromic. If it is it returns true.
	var normalExp = toMinimum(expressionString);
	var revcompExp = complementExp(reverseExp(normalExp));
	if (normalExp == revcompExp)	{
		return true;
	}
	else	{
		return false;
	}
}
function prepareAlign (alignArray)	{
	titleArray = new Array (alignArray.length - 1);
	sequenceArray = new Array (alignArray.length - 1);
	longestTitle = 0;
	for (var i = 1; i < alignArray.length; i++)	{
		titleArray[i-1] = alignArray[i].match(/[^\f\n\r]+[\f\n\r]/);
		titleArray[i-1] = filterFastaTitle(titleArray[i-1].toString());
		if (titleArray[i-1].length > longestTitle)	{
			longestTitle = titleArray[i-1].length;
		}
		sequenceArray[i-1] = alignArray[i].replace(/[^\f\n\r]+[\f\n\r]/,"");
		sequenceArray[i-1] = filterAlignSeq (sequenceArray[i-1]);
	}
	return true;
}
function randSeq (components,lengthOut)	{
	var sequence = "";
	var tempNum = 0;
	var tempChar = "";
	for (var j = 0; j < (lengthOut); j++)	{
		tempNum = (Math.random() * components.length);
		tempNum = Math.round(tempNum);
		if (tempNum == components.length)	{
			tempNum = 0;
		}
		tempChar = components[tempNum];
		sequence = sequence + tempChar;
		if (sequence.length == 60)	{
			outputWindow.document.write(sequence +'<br>\n');
			sequence = "";
		}
	}
	outputWindow.document.write(sequence + '<br>\n');
	return true;
}

//去掉fasta格式序列头的>行，只返回序列
function removeFastaTitleDna (dnaSequence)	{
	fastaSequenceTitle = "";
	if (dnaSequence.search(/\>[^\f\n\r]+[\f\n\r]/) != -1)	{
		fastaSequenceTitle = (dnaSequence.match(/\>[^\f\n\r]+[\f\n\r]/,"")).toString();
		fastaSequenceTitle = fastaSequenceTitle.replace(/\>|[\f\n\r]/g,"");
		fastaSequenceTitle = filterFastaTitle(fastaSequenceTitle);
		dnaSequence = dnaSequence.replace(/\>[^\f\n\r]+[\f\n\r]/,"");
	}
	return dnaSequence;
}
function removeFastaTitleProtein (proteinSequence)	{
	fastaSequenceTitle = "";
	if (proteinSequence.search(/\>[^\f\n\r]+[\f\n\r]/) != -1)	{
		fastaSequenceTitle = (proteinSequence.match(/\>[^\f\n\r]+[\f\n\r]/,"")).toString();
		fastaSequenceTitle = fastaSequenceTitle.replace(/\>|[\f\n\r]/g,"");
		fastaSequenceTitle = filterFastaTitle(fastaSequenceTitle);
		proteinSequence = proteinSequence.replace(/\>[^\f\n\r]+[\f\n\r]/,"");
	}
	return proteinSequence;
}
function realLength(aString)	{	//determines length of expression without square brackets and slashes.
	aString = aString.replace(/\[[a-z]+\]/gi,"X");
	aString = aString.replace(/\//g,"");
	return aString.length;
}
function replaceSpace (textToChange)	{
	return textToChange.replace(/ /g,"&nbsp;");
}
function reverse (dnaSequence)	{	
	var tempDnaArray = new Array;
	if (dnaSequence.search(/./) != -1)	{
		tempDnaArray = dnaSequence.match(/./g);
		tempDnaArray = tempDnaArray.reverse();
		dnaSequence = tempDnaArray.join("");
	}
	return dnaSequence;
}
function reverseExp (expressionString)	{	//reverses expression. Expression needs both slashes. Used in conjunction with complementExp function for non-palindromic expressions.
	var tempExp = "";
	var param = "";				//param holds the parameters located after the second slash. They are not reversed.
	var z = expressionString.length - 1;
	while (z > 0)	{
		if (expressionString.charAt(z) == "/")	{
			z = z - 1;
			while (expressionString.charAt(z) != "/")	{
				tempExp = tempExp + expressionString.charAt(z);
				z = z - 1;
			}
		}
		else	{
			param = expressionString.charAt(z) + param;
			z = z - 1;
		}
	}
	tempExp = tempExp.replace(/\[/g,"z");
	tempExp = tempExp.replace(/\]/g,"[");
	tempExp = tempExp.replace(/\z/g,"]");
	tempExp = "/" + tempExp + "/" + param;
	return tempExp;
}
function rightNum(theNumber,sequenceToAppend,lengthOfColumn,tabIn)	{
	var j = 0;
	var tempString = "";
	theNumber = theNumber.toString();
	for (var j = theNumber.length; j < lengthOfColumn; j++)	{
		tempString = tempString + "&nbsp;";
	}
	theNumber = tempString + theNumber + "&nbsp;";
	sequenceToAppend = sequenceToAppend + theNumber + tabIn;
	return sequenceToAppend;
}
function sequenceStats (sequence,arrayOfItems)	{	//arrayOFItems are regular expressions. A number included with each regular expression serves as a multiplier for the percentage calculation. Any additional text will appear next to the pattern when the results are given.
	var originalLength = sequence.length;
	outputWindow.document.write ('<table BORDER><tr><td ALIGN=LEFT>' + 'Pattern:' + '</td><td ALIGN=LEFT>' + 'Times Found:' + '</td><td ALIGN=LEFT>' + 'Percentage:' + '</td></tr>\n');
	for (var i = 0; i < (arrayOfItems.length); i++)	{
		var tempNumber = 0;
		var matchExp = arrayOfItems[i].match(/\/.+\//) + "gi";
		matchExp = eval(matchExp);
		if (sequence.search(matchExp) != -1)	{
			tempNumber = ((sequence.match(matchExp)).length);
		}
		outputWindow.document.write ('<tr><td ALIGN=LEFT>' + ((arrayOfItems[i].match(/\([^\(]+\)/)).toString()).replace(/\(|\)/g,"") + '</td><td ALIGN=LEFT>' + tempNumber + '</td><td ALIGN=LEFT>' + Math.round(100 * tempNumber * parseFloat(arrayOfItems[i].match(/\d+/))/originalLength) + '</td></tr>\n');
	}
	outputWindow.document.write('</table>\n');
	return true;
}
function shuffle (sequence)	{
	var tempSeq = "";
	var tempChar = "";
	var tempString1 = "";
	var tempString2 = "";
	var randNum = 0;
	var maxNum = 0;
	while (sequence.length > 0)	{
		maxNum = sequence.length;
		randNum = (Math.random() * maxNum);
		randNum = Math.round(randNum);
		if (randNum > (maxNum - 1)) {
			randNum = 0;
		}
		tempChar = sequence.charAt(randNum);
		tempSeq = tempSeq + tempChar;		
		tempString1 = sequence.substring(0,randNum);
		tempString2 = sequence.substring((randNum + 1),sequence.length);
		sequence = tempString1 + tempString2;
		if (tempSeq.length == 60)	{
			outputWindow.document.write(tempSeq +'<br>\n');
			tempSeq = "";
		}
	}
	outputWindow.document.write(tempSeq);
	return true;
}
function simpleStatsFasta (sequence)	{
	var stringToReturn = '<b>Results for ' + sequence.length + ' residue sequence ';
	if (fastaSequenceTitle.search(/[^\s]/) != -1)	{
		stringToReturn = stringToReturn + '"' + fastaSequenceTitle + '"';
	}
	stringToReturn = stringToReturn + ' starting "' + sequence.substring(0,10) + '".</b><br>\n';
	return stringToReturn;
}
//去掉dna序列[]中非gatcn的东西
function toMinimum(expressionString)	{//omits redundant information in a DNA regular expression.
	var z = 0;
	var dnaResidues = new Array ("g","a","t","c","n");
	var newString = "";
	while (z < expressionString.length)	{
		var foundArray = new Array ("","","","","");
		if (expressionString.charAt(z) == "[")	{
			while (expressionString.charAt(z) != "]")	{
				for (var i = 0; i < dnaResidues.length; i++)	{	
					if (expressionString.charAt(z).toUpperCase() == dnaResidues[i].toUpperCase())	{
						foundArray[i] = dnaResidues[i];
						break;
					}
				}
				z = z + 1;
			}
			newString = newString + "[";
			for (var i = 0; i < dnaResidues.length; i++)	{
				if (foundArray[i] != "")	{
					newString = newString + foundArray[i];
				}
			}
			newString = newString +  "]";
			z = z + 1;
		}
		else	{
			newString = newString + expressionString.charAt(z);
			z = z + 1;
		}
	}
	return newString;
}
//翻译
function translate (dnaSequence,startPos,strand,geneticCode)	{
	var proteinSequence = "";
	var codon = "";
	var proteinChar = "X";
	startPos = parseInt(startPos);
	if (strand == "reverse")	{
		dnaSequence = reverse(complement(dnaSequence));
	}
	for (var i = startPos; i < (dnaSequence.length - 2); i = i + 3)	{
		codon = dnaSequence.substring(i,(i+3));
		proteinChar = "X";
		for (var j = 0; j < geneticCodeMatchExp.length; j++)	{
			if (codon.search(geneticCodeMatchExp[j]) != -1)	{
				proteinChar = geneticCodeMatchResult[j];
				break;
			}
		}
		proteinSequence = proteinSequence + proteinChar;
	}
	return proteinSequence;
}
function translateUpperCase (mixedCaseDna)	{
	var stringToReturn = "";
	for (var n = 0; n < mixedCaseDna.length; n++)	{
		if (foundTransStart == false)	{
			if (mixedCaseDna.charAt(n).search(/[A-Z]/) != -1)	{
				foundTransStart = true;
				stringToReturn = stringToReturn + upperCaseTrans.charAt(0);
				upperCaseTrans = upperCaseTrans.substring(1,upperCaseTrans.length);
				aminoAcidsCounter = aminoAcidsCounter + 1;
				frameCounter = 4;
			}
			else	{
				stringToReturn = stringToReturn + " ";
			}
		}
		else	{
			if (mixedCaseDna.charAt(n).search(/[A-Z]/) == -1)	{
				stringToReturn = stringToReturn + " ";
			}
			else	{
				if (frameCounter % 3 != 0)	{
					frameCounter = frameCounter + 1;
					stringToReturn = stringToReturn + " ";
				}
				else	{
					stringToReturn = stringToReturn + upperCaseTrans.charAt(0);
					upperCaseTrans = upperCaseTrans.substring(1,upperCaseTrans.length);
					aminoAcidsCounter = aminoAcidsCounter + 1;
					frameCounter = frameCounter + 1;
				}
			}
		}
	}
	return stringToReturn;
}
function verifyDigits (theNumber)	{
	if (theNumber.search(/\d/) == -1) {
		alert ("Please enter a number");
		return false;
	}
}
//判断提交序列中是否有gatc外的另外字符，如有，则出提示。
function verifyDna (dnaSequence)	{
	if (dnaSequence.search(/[^gatcn\s]/i) != -1) {	
		alert ("The sequence contains non-DNA characters, which will be omitted.");
	}
	return true;
}
function verifyGenBank (genBankFile) {
	if ((genBankFile.search(/LOCUS/) == -1) || (genBankFile.search(/DEFINITION/) == -1) || (genBankFile.search(/ACCESSION/) == -1) || (genBankFile.search(/ORIGIN/) == -1)) {
		alert ("Please enter the contents of a GenBank file.");
		return false;
	}
	return true;
}
function verifyGenBankFeat (genBankFile) {
	if ((genBankFile.search(/LOCUS/) == -1) || (genBankFile.search(/DEFINITION/) == -1) || (genBankFile.search(/ACCESSION/) == -1) || (genBankFile.search(/ORIGIN/) == -1)) {
		alert ("Please enter the contents of a GenBank file.");
		return false;
	}
	if (genBankFile.search(/FEATURES {13}/) == -1)	{
		alert ("The file has no defined features.");
		return false;
	}
	return true;
}
function verifyProtein (proteinSequence)	{
	if (proteinSequence.search(/[^acdefghiklmnpqrstvwyz\*\s]/i) != -1) {	
		alert ("The sequence contains non-protein characters, which will be omitted.");
	}
	return true;
}