function validateForm() {
	const file = document.getElementById("file").files[0];
	const startSeq = document.getElementById("start_seq").value;
	const stopSeq = document.getElementById("stop_seq").value;
	const minCount = document.getElementById("min_count").value;
        if (!file || !startSeq || !stopSeq|| !minCount) {
		alert("Please fill out all fields.");
		return false;
	}

	if (!file.name.endsWith(".fastq")) {
		alert("Please upload a FASTQ file.");
		return false;
	}

	return true;
}
