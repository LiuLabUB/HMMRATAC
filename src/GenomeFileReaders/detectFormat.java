package GenomeFileReaders;

public class detectFormat {
	
	private String fileName;
	
	
	public detectFormat(String f){
		fileName = f.toLowerCase();
	}
	
	public boolean detectSAM(){
		if (fileName.endsWith("sam")){return true;}
		else{return false;}
	}
	public boolean detectBAM(){
		if (fileName.endsWith("bam")){return true;}
		else{return false;}
	}
	public boolean detectBED(){
		if (fileName.endsWith("bed")){return true;}
		else {return false;}
	}
	
	
}
