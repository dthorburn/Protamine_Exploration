suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))

gtf_path <- "/path/to/GTFs" ## Directory with all the GTFs downloaded from Ensembl metazoa. 
orths <- fread("FBgn0013300_Orthology.csv") ## Get from OrthoFinder

## A function to convert GTF coordinates into BED12 format and if there are no UTRs it creates synthetic ones by extending the UTRs by the mean ditance of the real ones in this data.
## Empirical estimates in this dataset are 5'=~200 & 3'=~300, rounded to the nearest 100. 
## To only extract BEDs with non-synthetic UTRs use syn_3/5_length = 0
Estimate_UTRs <- function(transcript_ID, Assembly_Accession = NA){
	message(paste("\tCalculating:", transcript_ID))
	gtf_trans <- subset(temp_gtf, transcript_id == transcript_ID)
	x3 <- subset(gtf_trans, type == "three_prime_utr")[,width] %>% sum
	x5 <- subset(gtf_trans, type == "five_prime_utr")[,width] %>% sum

	output <- data.table(transcript_ID = transcript_ID,
						Assembly_Accession = Assembly_Accession,
						x3_UTR_Length = ifelse(x3 == 0, NA, x3),
						x5_UTR_Length = ifelse(x5 == 0, NA, x5))

	return(output)
}

## A function to convert GTF coordinates into BED12 format and if there are no UTRs it creates synthetic ones by extending the UTRs by the mean ditance of the real ones in this data.
## Empirical estimates in this dataset are 5'=~200 & 3'=~300, rounded to the nearest 100. 
Make_BED12 <- function(transcript_ID, syn_5_length = 0, syn_3_length = 0, x5_only = FALSE, x3_only = TRUE){
	message(paste("\t\tExtracting:", transcript_ID))
	gtf_trans <- subset(temp_gtf, transcript_id == transcript_ID)
	jus_trans <- subset(gtf_trans, type == "transcript")
	jus_exons <- subset(gtf_trans, type == "exon")[,start]

	blocks <- grepl(pattern = "exon", gtf_trans$type)
	
	## Need to update to handle +/- strand
	utrs <- c(grepl(pattern = "five", gtf_trans$type) %>% sum >= 1,  grepl(pattern = "three", gtf_trans$type) %>% sum >= 1) %>% sum
	temp_name <- ifelse(utrs == 2, jus_trans$transcript_id, 
					ifelse(sum(grepl(pattern = "five", gtf_trans$type)) >= 1, paste0(jus_trans$transcript_id, "_syn3"),
						ifelse(sum(grepl(pattern = "three", gtf_trans$type)) >= 1, paste0(jus_trans$transcript_id, "_syn5"),
							ifelse(utrs == 0, paste0(jus_trans$transcript_id, "_synBoth"), warning(paste0(transcript_ID, ": error in UTR detection!"))))))

	bStarts_raw <- (jus_exons - jus_trans$start)
	bSizes_raw <- gtf_trans[blocks,width]
	start_raw <- (jus_trans$start - 1) ## It was not taking the first bp for some reason. This is what the converted USCS bed12's have instead. 
	end_raw <- jus_trans$end
	
	## Reversing the order if it's the negative strand - it's what USCS tools seem to do
	if(jus_trans$strand == "-"){
		bStarts_raw <- bStarts_raw %>% rev
		bSizes_raw <- bSizes_raw %>% rev
	}
	
	## Adding the synthetic UTRs if needed, and ensuring strandedness is taken into account. This got way more complicated than I had anticipated. Strandedness is a bitch. 
	## Turns out the same logic is applied as expected, but I needed to reverse the order of the blocks. 
	if(grepl(pattern = "syn5", temp_name)){
		message(paste0("\t\t\tAdding synthetic 5' to ", transcript_ID))
		if(jus_trans$strand == "+"){
			start_raw <- start_raw - syn_5_length
			bSizes_raw[1] <- bSizes_raw[1] + syn_5_length
			if(length(bStarts_raw) > 1){ 
				bStarts_raw[2:length(bStarts_raw)] <- bStarts_raw[2:length(bStarts_raw)] + syn_5_length
			}
		} else if(jus_trans$strand == "-"){
			end_raw <- end_raw + syn_5_length
			bSizes_raw[length(bSizes_raw)] <- bSizes_raw[length(bSizes_raw)] + syn_5_length
			if(length(bStarts_raw) > 1){ 
				bStarts_raw[1:(length(bStarts_raw)-1)] <- bStarts_raw[1:(length(bStarts_raw)-1)] + syn_5_length
			}
		}

	} 	else 	if(grepl(pattern = "syn3", temp_name)){
		message(paste0("\t\t\tAdding synthetic 3' to ", transcript_ID))
		if(jus_trans$strand == "+"){
			end_raw <- end_raw + syn_3_length
			bSizes_raw[length(bSizes_raw)] <- bSizes_raw[length(bSizes_raw)] + syn_3_length 
		} else if(jus_trans$strand == "-"){
			start_raw <- start_raw - syn_3_length
			bSizes_raw[1] <- bSizes_raw[1] + syn_3_length
			if(length(bStarts_raw) > 1){ 
				bStarts_raw[2:length(bStarts_raw)] <- bStarts_raw[2:length(bStarts_raw)] + syn_3_length
			}
		}

	}	else 	if(grepl(pattern = "Both", temp_name)){
		message(paste0("\t\t\tAdding synthetic 5' & 3' to ", transcript_ID))
		if(jus_trans$strand == "+"){
			start_raw <- start_raw - syn_5_length
			end_raw <- end_raw + syn_3_length
			bSizes_raw[1] <- bSizes_raw[1] + syn_5_length
			bSizes_raw[length(bSizes_raw)] <- bSizes_raw[length(bSizes_raw)] + syn_3_length
			if(length(bStarts_raw) > 1){ 
				bStarts_raw[2:length(bStarts_raw)] <- bStarts_raw[2:length(bStarts_raw)] + syn_5_length
			}
		} else if(jus_trans$strand == "-"){	
			start_raw <- start_raw - syn_3_length
			end_raw <- end_raw + syn_5_length
			bSizes_raw[1] <- bSizes_raw[1] + syn_3_length
			bSizes_raw[length(bSizes_raw)] <- bSizes_raw[length(bSizes_raw)] + syn_5_length
			if(length(bStarts_raw) > 1){ 
				bStarts_raw[2:length(bStarts_raw)] <- bStarts_raw[2:length(bStarts_raw)] + syn_3_length
			}
		}
	}

	temp_bed <- data.table(	chr = gtf_trans$seqnames %>% unique,
							start = start_raw,
							end = end_raw,
							name = temp_name, ## can be updated to include more information like the CDS start and end.
							score = 0,
							strand = jus_trans$strand,
							thickStart = start_raw,
							thickEnd = end_raw,
							itemRGB = 0,
							blockCount = blocks %>% sum,
							blockSizes =  bSizes_raw %>% paste(collapse = ",") %>% paste0(.,","),
							blockStarts = bStarts_raw %>% paste0(., collapse = ",") %>% paste0(.,","))
	## Removing negative values in the start column and annotating the name. Currently this does not handle 3'UTRs going beyond boundaries, but can if be if fais are added. Only 1 edge case here so I did it manually. 
	temp_bed[,"name" := ifelse(start <= 0, paste0(name, "_Adj_5UTR"), name)]
	temp_bed[,"start" := ifelse(start <= 0, 1, start)]
	temp_bed[,"thickStart" := ifelse(thickStart <= 0, 1, thickStart)]
	
	if(x5_only == FALSE & x3_only == FALSE){
		## Case 1: Not extracting just the UTRs
#		return(temp_bed)
	} 	else 	if(x5_only == TRUE){
		## Case 2: Extracting just the 5' UTRs 
		message(paste0("\t\t\t\tExtracting only 5' UTR bed for ", transcript_ID))
		## Case 2.1: Synthetic 5' 
		if(grepl(pattern = "Both|syn5", temp_name)){
			## Dealing with start and end locations			
			if(jus_trans$strand == "+"){
				## Case 2.1.1: Synthetic forward strand
				temp_bed$end <- jus_trans$start - 1  
			} else if(jus_trans$strand == "-"){
				## Case 2.1.2: Synthetic reverse strand
				temp_bed$start <- jus_trans$end
			}
			## Changing the number of bedBlocks to 1 since it's just the synthetic UTR to be extracted. 
			## and then altering the rest of the temp bed
			temp_bed$blockCount <- 1
			temp_bed$name <- paste0(temp_bed$name, "_5only")
			temp_bed$blockSizes <- c(syn_5_length, "") %>% paste(collapse = ",")
			temp_bed$blockStarts <- c(0, "") %>% paste(collapse = ",")
			temp_bed$thickStart <- temp_bed$start
			temp_bed$thickEnd <- temp_bed$end
		} else {
			## Taking the coordinates of the UTR from the GTF file and sorting it
			gtf_trans_utr <- subset(gtf_trans, type == "five_prime_utr")
			setkey(gtf_trans_utr, start, end)

			## Adjusting the temporary bed
			temp_bed$start <- gtf_trans_utr$start[1] - 1
			temp_bed$end <- gtf_trans_utr$end[nrow(gtf_trans_utr)]
			temp_bed$name <- paste0(temp_bed$name, "_5only")
			temp_bed$blockCount <- grepl(pattern = "five", gtf_trans$type) %>% sum
			temp_bed$blockSizes <- c(gtf_trans_utr$width, "") %>% paste(collapse = ",")
			temp_bed$blockStarts <- c((gtf_trans_utr$start - gtf_trans_utr$start[1]), "") %>% paste(collapse = ",")
			temp_bed$thickStart <- temp_bed$start
			temp_bed$thickEnd <- temp_bed$end
		}
#		return(temp_bed)
		## Logic here is to extract the first/last n number of entries in the blockSizes and blockStarts variables of the bed to create the UTR only entries. 
		## The difficulty will be to account for the strandedness and when syntehtic ones are included. 
		## Another problem is using the "exon" tag above, as some exons with UTRs attached (i.e., no gap combine exon coordinates as CDS+UTR; see TDAL004025-RA for example)
		## The above code DOES NOT ACCOUNT FOR THE STRANDEDNESS OR EXON PROBLEMS - FIX ME!

	}	 else 	if(x3_only == TRUE){
		## Case 3: extracting just the 3'
		message(paste0("\t\t\t\tExtracting only 3' UTR bed for ", transcript_ID))
		gtf_trans_utr <- subset(gtf_trans, type == "three_prime_utr")
		setkey(gtf_trans_utr, start, end)
		
		## Case 3.1: Synthetic 3'
		if(grepl(pattern = "Both|syn3", temp_name)){
			if(jus_trans$strand == "+"){
				## Case 3.1.1: Synthetic forward strand
				temp_bed$start <- jus_trans$end # + 1 - Is this necessary? 
			} else if(jus_trans$strand == "-"){
				## Case 3.1.2: Synthetic reverse strand
				temp_bed$end <- jus_trans$start - 1
			}
			temp_bed$name <- paste0(temp_bed$name, "_3only")
			temp_bed$blockCount <- 1
			temp_bed$blockSizes <- c(syn_3_length, "") %>% paste(collapse = ",")
			temp_bed$blockStarts <- c(0, "") %>% paste(collapse = ",")
			temp_bed$thickStart <- temp_bed$start
			temp_bed$thickEnd <- temp_bed$end
		} else {
			## Case 3.2: Real 3' - The same code as the real 5', but logically I'm going to keep it organised like this as it's easier to follow.
			gtf_trans_utr <- subset(gtf_trans, type == "three_prime_utr")
			setkey(gtf_trans_utr, start, end)
			temp_bed$start <- gtf_trans_utr$start[1] - 1
			temp_bed$end <- gtf_trans_utr$end[nrow(gtf_trans_utr)]
			temp_bed$name <- paste0(temp_bed$name, "_3only")
			old_blockCount <- temp_bed$blockCount
			temp_bed$blockCount <- grepl(pattern = "three", gtf_trans$type) %>% sum
			temp_bed$blockSizes <- c(gtf_trans_utr$width, "") %>% paste(collapse = ",")
			temp_bed$blockStarts <- c((gtf_trans_utr$start - gtf_trans_utr$start[1]), "") %>% paste(collapse = ",")
			temp_bed$thickStart <- temp_bed$start
			temp_bed$thickEnd <- temp_bed$end
		}
#		return(temp_bed)
	}
	## This is to handle the cases when synthetic primers are not present and UTR only options are TRUE. 
	## It will not return the empty beds
	if(grepl(pattern = "^0,$", temp_bed$blockSizes)){
		message(paste0("\t\t\t\tWARN: No UTRs to extract - not returning bed for: ", temp_name))
	} else {
		return(temp_bed)
	}
}
	
## A loop to go through the assembly IDs and take all 
for(spec in unique(orths$Assembly_Accession)){
	move_to_next <- 1
	message(paste("Starting on:", spec))
	
	## Subsetting the candidates
	temp_orths <- subset(orths, Assembly_Accession == spec)

	## Reading in species data
	temp_gtf_path <- grep(pattern = spec, dir(gtf_path), value = TRUE) %>% paste0(gtf_path, "/", .)
	tryCatch(
		temp_gtf <- rtracklayer::import(temp_gtf_path) %>% as.data.table, ## Annoying this didn't work on the same line for whatever reason 
		error=function(e){ move_to_next <<- 0; error_no <<- 1 }
	)

	if(move_to_next == 1){
		temp_utrs <- mapply(FUN = Estimate_UTRs, temp_orths$Transcript_ID, spec, USE.NAMES = FALSE) %>% t %>% as.data.table ## This format is problematic. Save and relad as DT to get around lists. 
		temp_bed <- lapply(temp_orths$Transcript_ID, FUN = Make_BED12) %>% rbindlist %>% unique
		fwrite(file = paste0(spec, "_protamine_orths_Synth.bed"), sep = "\t", temp_bed, col.names = FALSE)
		if(spec == unique(orths$Assembly_Accession)[1]){
			all_beds <- temp_bed
			all_utrs <- temp_utrs
		}	else 	{
			all_beds <- rbind(all_beds, temp_bed)
			all_utrs <- rbind(all_utrs, temp_utrs)
		}
	}	else 	if(move_to_next == 0){
		message(paste("WARN:", spec, "GTF failed to read in. Try running in alone to debug."))
	}

}
fwrite(file = "All_BEDS_wSynth.csv", all_beds)
fwrite(file = "All_UTRs_wSynth.csv", all_utrs)
