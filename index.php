<html>

<head>
	<title>CRISPR-A</title>
	<meta http-equiv="content-type" content="text/html; charset=utf-8" />
	<meta name="viewport" content="width=device-width, initial-scale=1.0">
	<meta name="description" content="" />
	<meta name="keywords" content="" />
	<link rel="preconnect" href="https://fonts.googleapis.com">
	<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
	<link href="https://fonts.googleapis.com/css2?family=Public+Sans:wght@400;700&display=swap" rel="stylesheet">
	<script src="js/jquery.min.js"></script>
	<script src="js/config.js"></script>
	<script src="js/skel.min.js"></script>
	<link rel="stylesheet" href="css/skel-noscript.css" />
	<link rel="stylesheet" href="css/normalise.css" />
	<link rel="stylesheet" href="css/style.css" />
	<link rel="stylesheet" href="css/style-desktop.css" />
	<link rel="stylesheet" href="css/style-tablet.css" />
	<link rel="stylesheet" href="css/style-mobile.css" />

	<!-- Global site tag (gtag.js) - Google Analytics -->
	<script async src="https://www.googletagmanager.com/gtag/js?id=G-4LW5V5M941"></script>
	<script>
  		window.dataLayer = window.dataLayer || [];
  		function gtag(){dataLayer.push(arguments);}
  		gtag('js', new Date());
  		gtag('config', 'G-4LW5V5M941');
	</script>

<script type="text/javascript">
//Initialating some variables
var i = 1;
var mockArray = [];

</script>

</head>

<body class="home">

		<style>
		/* DRAG & DROP MOCK */
	      .upload-container {
	        position: relative;
	      }
	      .upload-container input {
			border: 2px dashed #ccc;
			border-radius: 20px;
	        background: #f1f1f1;
	        outline: 2px dashed #1cf453;
	        outline-offset: -10px;
	        padding: 100px 0px 100px 250px;
	        text-align: center !important;
	        width: 100%;
	      }
	      .upload-container input:hover {
	        background: #ddd;
	      }
	      .upload-container:before {
	        position: absolute;
	        bottom: 20%;
	        left: 39%;
	        content: "Drop your mock files here";
	        color: #696969;
	        font-weight: 900;
	      }
	      .upload-btn {
	        margin-left: 25%;
	        padding: 7px 20px;
	      }
		 /* MOCK DROPDOWN */
		.drop {
		  /* position: relative;
		  display: flex; */
		  /* width: 120%; */
		  height: 50%;
		  border-radius: .25em;
		  overflow: hidden;
		  background-color: white;
		  border-color: white;
		}
/* INFO TOOLTIP */
.tooltip {
  position: relative;
  display: inline-block;
  border-bottom: 1px dotted black;
}
.tooltip .tooltiptext {
  visibility: hidden;
  width: 3000%;
  background-color: #555;
  color: #fff;
  text-align: center;
  border-radius: 6px;
  padding: 5px 0;
  position: absolute;
  z-index: 1;
  bottom: 125%;
  left: -1075%;
  margin-left: -60px;
  opacity: 0;
  transition: opacity 0.3s;
}
.tooltip .tooltiptext::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 50%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent transparent transparent;
}
.tooltip:hover .tooltiptext {
  visibility: visible;
  opacity: 1;
}
	    </style>

	<header class="header">
		<div class="header-logo">
			<img src="Images/crispr-a-logo.svg" alt="CRISPR-A">
		</div>
		<div class="header-nav">
			<a class="header-nav__link button button--black" href="documentation.html">Help</a>
			<a class="header-nav__link button button--white" href="">Reference</a>
		</div>
		<div class="header-menu">
			<div class="header-menu__icon">
				<span class="closed">+</span>
				<span class="opened">-</span>
			</div>
		</div>
	</header>

	<main>

		<section class="home-hero hero">

			<div class="hero__text">
				<h1 class="hero__headline headline type--black">CRISPR-A helps <br/>you from design <br/>to analysis. </h1>
				<h2 class="hero__subheadline subheadline type--black">Simulate gene editing outcomes in your target.</h2>
				<h2 class="hero__subheadline subheadline type--black">Assess the quality of gene editing experiments using next generation sequencing (NGS) data.</h2>
			</div>

		</section>
		
			<div class="form__block form__block--ngs">
		<div class="form__wrapper">
			<h3 class="form__subheadline subheadline">How should we analyse your reads?</h3>
			<div class="input__design input__design--start">
				<span class="input__inner form__label">Type of Analysis:</span>
				<div class="input__inner form__choices">
					<div class="form__choice">
						<input type="radio" id="singlendCheckbox" name="design" >
						<label for="singlendCheckbox" class="form__label form__choice">Single-end Reads</label>
					</div>
					<div class="form__choice">
						<input type="radio" id="pairendCheckbox" name="design" checked>
						<label for="pairendCheckbox" class="form__label form__choice">Paired-end Reads</label>
					</div>
				</div>
		</div>
	</div>
</div>

		<form id="form" class="form" name="inputForm" enctype="multipart/form-data" action="upload.php" method="post" onsubmit="return validateForm()" >
		
					<form class="form" name="inputFormMock" enctype="multipart/form-data"  method="post" >
					<div id="dynamic_field_2">
							<div class="form__block form__block--ngs">
									<div class="form__wrapper">
											<h3 class="form__subheadline subheadline"> Control noise adding a negative control <i><div class="tooltip"> <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-question-circle" viewBox="0 0 16 16">
											   <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14zm0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16z"/>
											   <path d="M5.255 5.786a.237.237 0 0 0 .241.247h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286zm1.557 5.763c0 .533.425.927 1.01.927.609 0 1.028-.394 1.028-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94z"/>
											 </svg>
											   <span class="tooltiptext">Optional</span>
											 </div></i></h3>
											 <label class="form__label form__label--full">Drop your R1 files <i><div class="tooltip"> <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-question-circle" viewBox="0 0 16 16">
											 <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14zm0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16z"/>
											 <path d="M5.255 5.786a.237.237 0 0 0 .241.247h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286zm1.557 5.763c0 .533.425.927 1.01.927.609 0 1.028-.394 1.028-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94z"/>
											 </svg>
											 <span class="tooltiptext">Press ctrl while selecting your files to upload several at the same time</span>
											 </div></i></label>
											<div class="upload-container">
											  <input id="mockfile" name="uploaded_mock[]" type="file" multiple />
											 <div id="mockList"></div>
											</div>
							<div class="form__block form__block--ngs">
										 <label class="form__label form__label--full">Drop your R2 files <i><div class="tooltip"> <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-question-circle" viewBox="0 0 16 16">
										 <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14zm0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16z"/>
										 <path d="M5.255 5.786a.237.237 0 0 0 .241.247h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286zm1.557 5.763c0 .533.425.927 1.01.927.609 0 1.028-.394 1.028-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94z"/>
										 </svg>
										 <span class="tooltiptext">Once uploaded, check that the each pair of files are printed in the same position. Visit Help for more info.</span>
										 </div></i></label>
											<div class="upload-container">
											  <input id="mockfile2" class="pair" name="uploaded_mock1[]" type="file" multiple />
											  <div id="mockListRev"></div>
											</div>
											<br>
										</div>
											<div id="errorMock">
										</div>
										</div>
        
			<div id="dynamic_field">

				<div class="form__block form__block--ngs">

					<div class="form__wrapper" id="row1">

						<h3 class="form__subheadline subheadline">Analyse NGS datasets</h3>

						<div class="input__wrapper">
							<div class="input__inner input__inner--half">
								<span class="form__label">Upload forward read:</span>
								<input id="file" accept=".fastq, .fastq.gz, .fastq.zip, .fastq.gzip, .fastq.tar" name="uploaded_file[]"  type="file" required />
							</div>

							<div class="input__inner input__inner--half">
								<span class="form__label">Upload reverse read (optional):</span>
								<input name="uploaded_file1[]" type="file" accept=".fastq, .fastq.gz, .fastq.zip, .fastq.tar, .fastq.gzip" />
                            </div>

							<div id="error"></div>

						</div>

						<div class="input__wrapper input__wrapper--start">

							<span class="input__inner form__label">Reference genome:</span>
							<div class="input__inner form__choices">
								<div class="form__choice">
									<input type="radio" name="genome[]" id="human" value="human" class="required">
									<label for="human" class="form__label form__choice">Human </label>
								</div>
								<div class="form__choice">
									<input type="radio" id="mouse" name="genome[]" value="mouse">
									<label for="mouse" class="form__label form__choice">Mouse </label>
								</div>
								<div class="form__choice">
									<input type="radio" id="other" name="genome[]" value="other" checked>
									<label for="other" class="form__label form__choice">Other </label>
								</div>
							</div>

						</div>

						<div class="input__wrapper">
							<textarea id="reference" class="textarea text required" name="original[]" maxlength="500000" rows="5" placeholder="Paste reference sequence to be edited"></textarea>
						</div>
						<div id="error_reference"></div>

						<div class="input__wrapper">
			                <br>
							<div class="input__inner input__inner--half">
								<label for="protospacer" class="form__label form__label--full">Protospacer sequence:</label>
								<input id="protospacer" class="input__text required" name="grna[]" type="text" maxlength="40" value = "NNNNNNNNNNNNNNNNNNNN" /><br>
							</div>

							<div class="input__inner input__inner--half">

								<label class="form__label form__label--full">Target sequence (optional):</label>
								<input id="target" class="input__text" name="target[]" type="text" maxlength="100000"/><br>
                            	<br>
							</div>
                            
                            <div id="error_protospacer"></div>

						</div>

					<div class="input__wrapper">
						<div class="input__inner input__inner--half">
							<label class="form__label form__label--full">Cut site location in protospacer:</label>
							<input id="cutpos" class="input__text input__text required" name="cutsitePos[]" type="number" value="-3" /><br> 
						</div>

						<div class="input__inner input__inner--half">
							<label class="form__label form__label--full">Does this sample has a control mock sample?</label>
							<select style="width: 420px" id="mockRef" class="drop" name="mock_reference[]">
							  <option selected="">No</option>

								</select>
						</div>
						<div id="error_cutpos"></div>
					</div>
					</div>

					<div class="button__wrapper input__wrapper input__wrapper--end">
						<button type="button" name="add" id="add" class="cta button button--black">Add More</button>
					</div>

				</div>

			   </div>

			<div class="form__block form__block--presice">

				<div class="form__wrapper">

					<h3 class="form__subheadline subheadline">Do you need a highly precise analysis?</h3>

					<div class="input__wrapper">

						<div class="input__inner input__inner--half flex">
							<span class="form__label">Uni-molecular barcodes:</span>
							<div class="form__choices form__choices--no-margin flex">
								<div class="form__choice">
									<input type="radio" name="umis" id="umisyes" value="true">
									<label for="umisyes" class="form__label form__choice flex">Yes </label>
								</div>
								<div class="form__choice">
									<input type="radio" name="umis" id="umisno" value="false" checked>
									<label for="umisno" class="form__label form__choice flex">No </label>
								</div>
							</div>
						</div>

						<div class="input__inner input__inner--half flex">
							<span class="form__label">Spikes:</span>
							<div class="form__choices form__choices--no-margin flex">
								<div class="form__choice">
									<input type="radio" name="spike" id="spikeyes" value="yes">
									<label for="spikeyes" class="form__label form__choice flex">Yes </label>
								</div>
								<div class="form__choice">
									<input type="radio" name="spike" id="spikeno" value="no" checked>
									<label for="spikeno" class="form__label form__choice flex">No </label>
								</div>
							</div>
						</div>
					</div>

				</div>
				<div class="button__wrapper">
					<div id="errorMockRef">
					</div>
					<br>
					<br>
					<div id="errorGeneral">
					</div>
		
						<button class="cta button button--full button--green" type="submit" name="submit" value="Upload">Analyse</button>
					</div>
				</div>

				<div class="button__wrapper">
					<a style="text-align: center;" class="cta button button--full button--green" href="https://synbio.upf.edu/crispr-a/RUNS/tmp_482360152/">Example results</a>
				</div>


			</div>

			<div class="form__block form__block--grey">

				<div class="form__wrapper form__wrapper--white">
					<h3 class="form__subheadline subheadline">Analyse simulated dataset</h3>
					<textarea class="textarea textarea--grey text" name="simgeseq" maxlength="500000" rows="4" placeholder="Paste reference sequence to be edited"></textarea>

					<div class="input__wrapper">
						<div class="input__inner input__inner--half">
							<label class="form__label form__label--full">Cut site location:</label>
							<input class="input__text input__text--grey" name="cutsite" type="number"/><br>
						</div>
					</div>

					<div class="input__wrapper">
						<div class="input__inner input__inner--half">
							<label class="form__label form__label--full">Protospacer sequence:</label>
							<input class="input__text input__text--grey" name="grnaSim" type="text" maxlength="40" value = "NNNNNNNNNNNNNNNNNNNN"/><br>
						</div>
					</div>

				</div>

				<div class="button__wrapper button__wrapper--last input__wrapper input__wrapper--end ">
					<button type="submit" name="submit" value="Upload" class="cta button button--black">Get simulation</button>
				</div>

			</div>

		</form>

				<!-- <button onclick="location.href='documentation.html'">Help</button> -->



		<!-- <h4>Reference: </h4>
		<a href="http://bioinformatics.oxfordjournals.org/content/early/2014/06/30/bioinformatics.btu427.full.pdf">Guell M, Yang L, Church G (2014) Genome Editing Assessment using CRISPR Genome Analyzer (CRISPR-GA) <i>Bioinformatics.</i></a> -->

	</main>

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
            <p class="footer__text text type--white">Marta Sanvicente Garc�a:</p>
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
                    <a href="TermsOfService.pdf" class="footer-copyright__text footer__link type--dark-grey">Legal Notice</a>
                    <a href="PrivacyPolicy.pdf" class="footer-copyright__text footer__link footer__link--last type--dark-grey">Privacy Policy</a>
                </li>
            </ul>
        </div>

        <div class="footer__column">
            <ul class="footer__list">
                <li class="footer__item">
                    <span class="footer-copyright__text footer-copyright__link type--dark-grey">CRISPR-A 2021</span>
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


	<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
	<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.2/js/bootstrap.min.js"></script>

	</body>
	</html>



<script type="text/javascript">
//SINGLE END PAIRED END BLOCK SWITCH
$('#pairendCheckbox').on('click',function() {
	var cb = $('#pairendCheckbox').is(':checked');
	$('.pair').prop('disabled', !cb);
});
$('#singlendCheckbox').on('click',function() {
	var cb = $('#singlendCheckbox').is(':checked');
	$('.pair').prop('disabled', cb);
	errorMock.innerText = "";
	mockListRev.innerText = "";
});

const protospacer = document.getElementById('protospacer')
const reference = document.getElementById('reference') 
const target = document.getElementById('target')
const file = document.getElementById('file')
const form = document.getElementById('form')
const cutpos = document.getElementById('cutpos') 
 
//error pop up's for Reads 
const errorElement = document.getElementById('error')

//error pop up's for Reference
const errorRef = document.getElementById("error_reference")

//error pop up's for protospacer
const errorProto = document.getElementById("error_protospacer")

//error pop up's for cut site position
const errorCut = document.getElementById("error_cutpos") 

form.addEventListener('submit', (e) => {
 	errorMockRef.innerText = ""
	errorElement.innerText = ""
    errorRef.innerText = ""
    errorProto.innerText = ""
	errorGeneral.innerText = ""
  	errorCut.innerText = ""
  //PROTOSPACEER CONDITIONS
  let messages_prot = [] 
  if (protospacer.value === '' || protospacer.value === null) {
    messages_prot.push('You must provide a protospacer sequence')
  }else{
    if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(protospacer.value)!==true) {
     messages_prot.push("The protospacer doesn't look like DNA to me")
    }
    if (protospacer.value.length < 10 || protospacer.value.length > 30) {
        messages_prot.push("Are you sure this is the protospacer?")
    }
  }

  //CUT SITE CONDITIONS
  let message_cut = []  //modif
  if(Math.abs(cutpos.value) > protospacer.value.length){  
	  message_cut.push("Please select a value within the length of the protospacer") 
  }

  //TARGET CONDITIONS
  if (target.value !== ''){
	if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(target.value)!==true) {
		messages_prot.push("The Target sequence doesn't look like DNA to me")
	}
	if (target.value.length < 20) {
			messages_prot.push("Target too small minimum of 20nt is required")
		}
	if (target.value.length > 500) {
		messages_prot.push("Target too big maximum of 500nt is accepted")
	}
  }
  
  //REFERENCE CONDITIONS
  const choice = document.getElementById('other').checked
  
  let messages_ref = []

  if (reference.value === '' || reference.value == null) {
	  if (choice===true){
		messages_ref.push('You must provide a refrence sequence')
	  }
  }else{
	
	if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(reference.value)!==true) {
	messages_ref.push("Doesn't look like DNA to me")
	}
	if (reference.value.length < 100) {
		messages_ref.push("Too short for the refernce minimum of 100nt is required")
	}
	if (reference.value.length > 5000) {
		messages_ref.push("Too big for the refernce maximum of 5000nt is accepted")
	}
	
	}

  //READS FILES
  let messages = []

  if (file.value === '' | file.value == null ) {
    messages.push("Having the forward read would be nice :)")
  }else{
	const fi = document.getElementById('file');
  	const fsize = fi.files.item(0).size;
	const size = Math.round((fsize /(1024*1024)));
	
	if (size >= 2000) {
    	messages.push("File too Big, please select a file less than 2Gb");
    } else if (size < 0.000002 ) {
        messages.push("File too small, please select a file greater than 2bytes");
        }
  }
  
    // MOCK USE CEHCK
  console.log(mockFiles);
  var difference = mockFiles.filter(x => mockReferences.indexOf(x) === -1);
  console.log(difference);
  if (difference.length !== 0) {
	 errorMockRef.innerText = "";
	 messages_mock_ref.push("The following mock files are not referenced to any treated sample: ")
	 console.log("The following mock files are not referenced to any treated sample:");
	 for(let h = 0; h < difference.length; h++){
		messages_mock_ref.push(difference[h] + ' ')
		console.log(difference[h]);
	 }
  }

  
 if (messages_prot.length > 0 || messages_ref.length > 0 || messages.length > 0 || messages_mock_ref.length > 0 || messages_mock_ref.length > 0 || message_cut.length > 0) {
	 e.preventDefault()
	 errorElement.innerText = messages.join(', ')
	 errorElement.style.color = '#ff632c';
	 
	 errorRef.innerText = messages_ref.join(', ')
	 errorRef.style.color = '#ff632c';
	 
	 errorProto.innerText = messages_prot.join(', ')
	 errorProto.style.color = '#ff632c';
	 
	 errorMockRef.innerText = messages.join(', ')
	 errorMockRef.style.color = '#ff632c';
	 
	 errorGeneral.innerText = "Ups, something seems to not add up, scroll up to double check"
	 errorGeneral.style.color = '#ff632c';
	 
	errorCut.innerText = message_cut.join(', ')  
	errorCut.style.color = '#ff632c';  
  }
})

var mockNum = 0;
var mockRevNum = 0;
var mockFilesRev=[];
var mockFiles=[];
var errorMock = document.getElementById("errorMock");
var mockList = document.getElementById("mockList");
$("#mockfile").on("change", function(){
  let messages_mock = [];
  let mockFiles = [];
  mockList.innerText = "";
  var f = document.getElementById('mockfile');
  var g = document.getElementById('mockfile2');
  var txt = "";
  console.log(f.value);
  $('.drop').empty();
  txt = "<option class='file'>No</option>";
  $('.drop').append(txt);
  if ('files' in f){
	for (var i = 0; i < f.files.length; i++) {
		var file = f.files[i];
		mockFiles.push(file.name);
		if ('name' in file) {
			txt = "<option class='file'>" + file.name + "</option>";
			$('.drop').append(txt);
		}
	}
	mockFiles.sort();
	var mockNum = f.files.length;
	var mockRevNum = g.files.length;
	for (var i = 0; i < mockNum; i++) {
		mockList.innerText = mockFiles.join(', ');
		mockList.style.fontStyle = "italic";
	}
	mockList.innerHTML = mockFiles;
	var checkBox = document.getElementById("pairendCheckbox");
	if (checkBox.checked == true){
		if (mockRevNum != mockNum) {
  		  messages_mock.push("Did you forget any file? The number of forward reads files does not match the reverse ones");
  		  } else {
 			 errorMock.innerText = "";
 		}
  	}
	if (messages_mock.length > 0) {
		  errorMock.innerText = messages_mock.join(', ');
		  errorMock.style.color = '#ff632c';
	  }
	}
  //Then the file list is cleaned to further repeat iterations and the options are added to the dropdown through a loop
  // $(".file").remove();
  // $('.drop').empty(); // Clean the drop-down
  // mockArray[0]= "<option class='file'>No</option>";
  // mockArray[1]=txt;
  // var a = mockArray.length;
  // for (var count = 0; count < a; count++) {
	// $('.drop').append(mockArray[count]);
	// }
});
$("#mockfile2").on("change", function(){
  let messages_mock = [];
  g = document.getElementById('mockfile2');
  f = document.getElementById('mockfile');
  // console.log(g.value);
  var mockNum = f.files.length;
  var mockRevNum = g.files.length;
	let mockFilesRev = [];
  if ('files' in g){
	  for (var i = 0; i < g.files.length; i++) {
		  var file = g.files[i];
		  mockFilesRev.push(file.name);
	  }
	  mockFilesRev.sort();
	  // const mockRevNum = mockFilesRev.length;
	  for (var i = 0; i < mockRevNum; i++) {
		  mockListRev.innerText = mockFilesRev.join(', ');
		  mockListRev.style.fontStyle = "italic";
	  }
	  if (mockRevNum != mockNum) {
		  messages_mock.push("Did you forget any file? The number of forward reads files does not match the reverse ones");
	  } else {
		  errorMock.innerText = "";
	  }
  }
  if (messages_mock.length > 0) {
	 // e.preventDefault();
	 errorMock.innerText = messages_mock;
	 errorMock.style.color = '#ff632c';
  }
  });
$(document).ready(function(){
  var i = 1;
  var mockReferences=[];
  var mockFiles=["No"];
  $("#add").on("click", function(){
	 i++;
	 	 $('#dynamic_field').append('<div class="form__block form__block--ngs"><div class="form__wrapper" id="row'+i+'"><h3 class="form__subheadline subheadline">Analyse NGS datasets</h3><div class="input__wrapper"><div class="input__inner input__inner--half"><span class="form__label">Upload forward read:</span><input id="file_'+i+'" type="file" name="uploaded_file[]" type="file" accept=".fastq, .fastq.gz, .fastq.zip, .fastq.tar" required /></div>  <div class="input__inner input__inner--half"><span class="form__label">Upload reverse read:</span><input name="uploaded_file1[]" class="pair" type="file" accept=".fastq, .fastq.gz, .fastq.zip, .fastq.tar" /></div> <div id="error_'+i+'"></div> </div><div class="input__wrapper input__wrapper--start"><span class="input__inner form__label">Reference genome:</span><div class="input__inner form__choices"><div class="form__choice"><input type="radio" name="genome[]'+i+'" id="human'+i+'" value="human"><label for="human'+i+'" class="form__label form__choice">Human </label></div><div class="form__choice"><input type="radio" id="mouse'+i+'" name="genome[]'+i+'" value="mouse"><label for="mouse'+i+'" class="form__label form__choice">Mouse </label></div><div class="form__choice"><input type="radio" id="other'+i+'" name="genome[]'+i+'" value="other" checked><label for="other'+i+'" class="form__label form__choice">Other </label></div></div></div><div class="input__wrapper"><textarea id="reference_'+i+'" class="textarea text" name="original[]" maxlength="500000" rows="5" placeholder="Paste reference sequence to be edited"></textarea></div> <div id="error_reference_'+i+'"></div> <div class="input__wrapper"><div class="input__inner input__inner--half"><label for = "protospacer" class="form__label form__label--full">Protospacer sequence:</label><input id="protospacer_'+i+'" class="input__text" name="grna[]" type="text" maxlength="40" value = "NNNNNNNNNNNNNNNNNNNN"/><br></div> <div class="input__inner input__inner--half"><label class="form__label form__label--full">Target sequence (optional):</label><input id="target_'+i+'" class="input__text" name="target[]" type="text" maxlength="100000"/><br> <br></div>  <div id="error_protospacer_'+i+'"> </div>  </div> <div class="input__wrapper"><div class="input__inner input__inner--half"> <label class="form__label form__label--full">Cut site location in protospacer:</label> <input id="cutpos'+i+'" class="input__text input__text required" name="cutsitePos[]'+i+'" type="number" value="-3" /><br> </div> <div class="input__inner input__inner--half"><label class="form__label form__label--full">Does this sample has a control mock sample?</label><select style="width: 420px" id="mockRef_'+i+'" class="drop" name="mock_reference[]'+i+'"><option selected="">No</option></select></div></div></div><div class="button__wrapper input__wrapper input__wrapper--end"><button type="button" name="remove" id="'+i+'" class="button button--white btn_remove">Remove</button><button type="button" name="add" id="add" class="add button button--black">Add More</button></div></div>')


		 // Here we capture the uploaded files and we print them as options

	 const protospacer = document.getElementById('protospacer_'+i)
	 const errorProto = document.getElementById("error_protospacer_"+i)
	 const target = document.getElementById('target_'+i)
	 const reference = document.getElementById("reference_"+i)
	 const errorRef = document.getElementById("error_reference_"+i)
	 const errorElement = document.getElementById("error_"+i)
	 const cutpos = document.getElementById('cutpos'+i) //modif
	 const errorCut = document.getElementById("error_cutpos_"+i)
	 var select = document.getElementById("mockRef_"+i);

	 var text = select.options[select.selectedIndex].text;
	 mockReferences.push(text);
	 let messages_mock_ref = []

	 form.addEventListener('submit', (e) => {
	 errorElement.innerText = ""
	 errorGeneral.innerText = ""
	 errorRef.innerText = ""
	 errorProto.innerText = ""
	 errorMockRef.innerText = ""
	 //PROTOSPACER CONDITIONS
	 let messages_prot = []
     if (protospacer.value === '' || protospacer.value === null) {
       messages_prot.push('You must provide a protospacer sequence')
     }else{
       if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(protospacer.value)!==true) {
        messages_prot.push("The protospacer doesn't look like DNA to me")
       }
       if (protospacer.value.length < 10 || protospacer.value.length > 30) {
           messages_prot.push("Are you sure this is the protospacer?")
       }
     }

	 //CUT SITE CONDITIONS	
	let message_cut = []  //modif	
	if(Math.abs(cutpos.value) > protospacer.value.length){  	
	   message_cut.push("Please select a value within the length of the protospacer") 	
	}
	
     //TARGET CONDITIONS
     if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(target.value)!==true && target.value !== '') {
       messages_prot.push("The Target sequence doesn't look like DNA to me")
     }

    //REFERENCE CONDITIONS
	const choice = document.getElementById('other'+i).checked
	 let messages_ref = []
     if (reference.value === '' || reference.value == null) {
		if (choice===true){
		messages_ref.push('You must provide a refrence sequence')
	  }
     }else{
       if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(reference.value)!==true) {
        messages_ref.push("Doesn't look like DNA to me")
       }
       if (reference.value.length < 100) {
           messages_ref.push("Too short for the refernce minimum of 100nt is required")
       }
     }

     //READS FILES
     let messages = []
	 if (file.value === '' | file.value == null ) {
		messages.push('Having the forward read would be nice :)')
	 }else{
	 const fi = document.getElementById('file');
  	 const fsize = fi.files.item(0).size;
	 const size = Math.round((fsize /(1024*1024)));

	if (size >= 2000) {
    	messages.push("File too Big, please select a file less than 2Gb");
    } else if (size < 0.000002 ) {
        messages.push("File too small, please select a file greater than 2mb");
        }
    }

	// MOCK USE CHECK
	console.log(mockFiles);
	var difference = mockFiles.filter(x => mockReferences.indexOf(x) === -1);
	console.log(difference);
	if (difference.length !== 0) {
	   errorMockRef.innerText = "";
	   messages_mock_ref.push("The following mock files are not referenced to any treated sample: ")
	   console.log("The following mock files are not referenced to any treated sample:");
	   for(let h = 0; h < difference.length; h++){
		   messages_mock_ref.push(difference[h] + ' ')
		   console.log(difference[h]);
	   }
	}

	if (messages_prot.length > 0 || messages_ref.length > 0 || messages.length > 0 || messages_mock_ref.length > 0) {
	e.preventDefault()
	errorElement.innerText = messages.join(', ')
	errorElement.style.color = '#ff632c';
	errorRef.innerText = messages_ref.join(', ')
	errorRef.style.color = '#ff632c';
	errorProto.innerText = messages_prot.join(', ')
	errorProto.style.color = '#ff632c';
	errorMockRef.innerText = messages.join(', ')
	errorMockRef.style.color = '#ff632c';
	errorCut.innerText = message_cut.join(', ')  	
	errorCut.style.color = '#ff632c';
	errorGeneral.innerText = "Ups, something seems to not add up, scroll up to double check"
	errorGeneral.style.color = '#ff632c';
  }
	 })
	 var f = document.getElementById('mockfile');
	 var txt = [];
	 console.log(mockReferences.length);
		 if ('files' in f){
		   for (var i = 0; i < f.files.length; i++) {
			   var file = f.files[i];
			   if ('name' in file) {
				   mockFiles.push(file.name);
				   txt.push("<option class='file'>" + file.name + "</option>");
			   }
		   }
		   let txtLen = txt.length;
		   for (let i = 0; i < txtLen; i++) {
			   $('.drop').append(txt[i]);
			 }
		 }

	 var cb = $('#singlendCheckbox').is(':checked');
	 $('.pair').prop('disabled', cb);
 });

 $(document).on("click",".add", function () {
   i++;
	 $('#dynamic_field').append('<div class="form__block form__block--ngs"><div class="form__wrapper" id="row'+i+'"><h3 class="form__subheadline subheadline">Analyse NGS datasets</h3><div class="input__wrapper"><div class="input__inner input__inner--half"><span class="form__label">Upload forward read:</span><input id="file_'+i+'" type="file" name="uploaded_file[]" type="file" accept=".fastq, .fastq.gz, .fastq.zip, .fastq.tar" required /></div>  <div class="input__inner input__inner--half"><span class="form__label">Upload reverse read:</span><input name="uploaded_file1[]" class="pair" type="file" accept=".fastq, .fastq.gz, .fastq.zip, .fastq.tar" /></div> <div id="error_'+i+'"></div> </div><div class="input__wrapper input__wrapper--start"><span class="input__inner form__label">Reference genome:</span><div class="input__inner form__choices"><div class="form__choice"><input type="radio" name="genome[]'+i+'" id="human'+i+'" value="human"><label for="human'+i+'" class="form__label form__choice">Human </label></div><div class="form__choice"><input type="radio" id="mouse'+i+'" name="genome[]'+i+'" value="mouse"><label for="mouse'+i+'" class="form__label form__choice">Mouse </label></div><div class="form__choice"><input type="radio" id="other'+i+'" name="genome[]'+i+'" value="other" checked><label for="other'+i+'" class="form__label form__choice">Other </label></div></div></div><div class="input__wrapper"><textarea id="reference_'+i+'" class="textarea text" name="original[]" maxlength="500000" rows="5" placeholder="Paste reference sequence to be edited"></textarea></div> <div id="error_reference_'+i+'"></div> <div class="input__wrapper"><div class="input__inner input__inner--half"><label for = "protospacer" class="form__label form__label--full">Protospacer sequence:</label><input id="protospacer_'+i+'" class="input__text" name="grna[]" type="text" maxlength="40" value = "NNNNNNNNNNNNNNNNNNNN"/><br></div> <div class="input__inner input__inner--half"><label class="form__label form__label--full">Target sequence (optional):</label><input id="target_'+i+'" class="input__text" name="target[]" type="text" maxlength="100000"/><br> <br></div>  <div id="error_protospacer_'+i+'"> </div>  </div> <div class="input__wrapper"><div class="input__inner input__inner--half"> <label class="form__label form__label--full">Cut site location in protospacer:</label> <input id="cutpos'+i+'" class="input__text input__text required" name="cutsitePos[]'+i+'" type="number" value="-3" /><br> </div> <div class="input__inner input__inner--half"><label class="form__label form__label--full">Does this sample has a control mock sample?</label><select style="width: 420px" id="mockRef_'+i+'" class="drop" name="mock_reference[]'+i+'"><option selected="">No</option></select></div><div id="error_cutpos_'+i+'"></div></div><div id="error_cutpos_'+i+'"></div></div></div><div class="button__wrapper input__wrapper input__wrapper--end"><button type="button" name="remove" id="'+i+'" class="button button--white btn_remove">Remove</button><button type="button" name="add" id="add" class="add button button--black">Add More</button></div></div>')
	 const protospacer = document.getElementById('protospacer_'+i)
	const errorProto = document.getElementById("error_protospacer_"+i)
	const target = document.getElementById('target_'+i)
	const reference = document.getElementById("reference_"+i)
	const errorRef = document.getElementById("error_reference_"+i)
	const errorElement = document.getElementById("error_"+i)
	const cutpos = document.getElementById('cutpos'+i) //modif
	const errorCut = document.getElementById("error_cutpos_"+i)
	var select = document.getElementById("mockRef_"+i);
	var text = select.options[select.selectedIndex].text;
	
	mockReferences.push(text);

	form.addEventListener('submit', (e) => {
	errorElement.innerText = ""
	errorRef.innerText = ""
	errorProto.innerText = ""
	errorGeneral.innerText = ""

	//PROTOSPACEER CONDITIONS
	let messages_prot = []
	if (protospacer.value === '' || protospacer.value === null) {
	  messages_prot.push('You must provide a protospacer sequence')
	}else{
	  if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(protospacer.value)!==true) {
	   messages_prot.push("The protospacer doesn't look like DNA to me")
	  }
	  if (protospacer.value.length < 10 || protospacer.value.length > 30) {
		  messages_prot.push("Are you sure this is the protospacer?")
	  }
	}

	//CUT SITE CONDITIONS
   	 let message_cut = []
     	 if(Math.abs(cutpos.value) > protospacer.value.length){
   		message_cut.push("Please select a value within the length of the protospacer")
    	 }

	//TARGET CONDITIONS 
	if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(target.value)!==true && target.value !== '') {
	  messages_prot.push("The Target sequence doesn't look like DNA to me")
	}

   //REFERENCE CONDITIONS
   const choice = document.getElementById('other'+i).checked
	let messages_ref = []
	if (reference.value === '' || reference.value == null) {
	if (choice===true){
	messages_ref.push('You must provide a refrence sequence')
  }
 }else{
   if(RegExp("^[A\A,C\C,T\T,N\N,G\G,a\a,c\c,g\g,t\t,n\n,' '\' ']+$").test(reference.value)!==true) {
	messages_ref.push("Doesn't look like DNA to me")
   }
   if (reference.value.length < 100) {
	   messages_ref.push("Too short for the refernce minimum of 100nt is required")
   }
 }


 //READS FILES
 let messages = []
 if (file.value === '' | file.value == null ) {
	messages.push('Having the forward read would be nice :)')
 }else{
 const fi = document.getElementById('file');
 const fsize = fi.files.item(0).size;
 const size = Math.round((fsize /(1024*1024)));

if (size >= 2000) {
	messages.push("File too Big, please select a file less than 2Gb");
} else if (size < 0.000002 ) {
	messages.push("File too small, please select a file greater than 2mb");
	}
}

// MOCK USE CEHCK
console.log(mockFiles);
var difference = mockFiles.filter(x => mockReferences.indexOf(x) === -1);
console.log(difference);
if (difference.length !== 0) {
	errorMockRef.innerText = "";
   messages_mock_ref.push("The following mock files are not referenced to any treated sample: ")
   console.log("The following mock files are not referenced to any treated sample:");
   for(let h = 0; h < difference.length; h++){
	   messages_mock_ref.push(difference[h] + ' ')
	   console.log(difference[h]);
   }
}

if (messages_prot.length > 0 || messages_ref.length > 0 || messages.length > 0 || messages_mock_ref.length > 0  || mesage_cut.length > 0) {
e.preventDefault()
errorElement.innerText = messages.join(', ')
errorElement.style.color = '#ff632c';
errorRef.innerText = messages_ref.join(', ')
errorRef.style.color = '#ff632c';
errorProto.innerText = messages_prot.join(', ')
errorProto.style.color = '#ff632c';
errorMockRef.innerText = messages.join(', ')
errorMockRef.style.color = '#ff632c';
errorCut.innerText = message_cut.join(', ')
errorCut.style.color = '#ff632c';
errorGeneral.innerText = "Ups, something seems to not add up, scroll up to double check"
errorGeneral.style.color = '#ff632c';
}

 })
	var f = document.getElementById('mockfile');
	var txt = [];
	if ('files' in f){
	  for (var i = 0; i < f.files.length; i++) {
		  var file = f.files[i];
		  if ('name' in file) {
			  txt.push("<option class='file'>" + file.name + "</option>");
		  }
	  }
	  let txtLen = txt.length;
	  for (let i = 0; i < txtLen; i++) {
		  $('.drop').append(txt[i]);
		}
	}
	var cb = $('#singlendCheckbox').is(':checked');
	$('.pair').prop('disabled', cb);
 	});



 $(document).on('click', '.btn_remove', function(){
   var button_id = $(this).attr("id");
 $('#row'+button_id+'').siblings('.button__wrapper').remove();
   $('#row'+button_id+'').remove();
 });

 $('.header-menu__icon').click(function(){
 	$('.header-nav').toggleClass('visible');
 	$(this).toggleClass('opened');
 });

 $("#submit").on('click',function(){
   var formdata = $("#add_name").serialize();
   $.ajax({
 	//url   :"action.php",
 	url : postURL,
 	method:"POST",  //added
 	//type  :"POST",
 	type:'json',
 	data  :formdata,
 	cache :false,
 	success:function(data){ //instead of result --> data
 	  alert(data);
 	  $("#add_name")[0].reset();
 	}
   });
 });
 });
 </script>