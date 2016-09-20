package HMMR_ATAC;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import Node.TagNode;

public class GenomeFileReader {

	private ArrayList<TagNode> _genMap;
	
	public GenomeFileReader(File input) throws FileNotFoundException{
		_genMap = new ArrayList<TagNode>();
		TagNode temp = null;
		Scanner inFile =new Scanner ((Readable) new FileReader(input));
		
		while(inFile.hasNext()){
			String chr = inFile.next();
			int size = inFile.nextInt();
			temp = new TagNode(chr,0,size);//Note the 20000 value is so it can be used in Bin.java
			_genMap.add(temp);
			
			
		}
		
	}
	public GenomeFileReader(String input) throws FileNotFoundException{
		_genMap = new ArrayList<TagNode>();
		TagNode temp = null;
		Scanner inFile =new Scanner ((Readable) new FileReader(input));
		
		while(inFile.hasNext()){
			String chr = inFile.next();
			int size = inFile.nextInt();
			temp = new TagNode(chr,0,size);//Note the 20000 value is so it can be used in Bin.java
			_genMap.add(temp);
			
			
		}
		
	}
	
	public ArrayList<TagNode> getMap(){
		return _genMap;
	}
	
	public ArrayList<TagNode> getFilteredMap(){
		ArrayList<TagNode> filteredMap = new ArrayList<TagNode>();
		TagNode temp = null;
		for (int i = 0;i < _genMap.size();i++){
			String chr = _genMap.get(i).getChrom();
			int start = _genMap.get(i).getStart();
			int stop = _genMap.get(i).getStop();
			int newStart = start + 20000;
			int newStop = stop - 20000;
			temp = new TagNode(chr,newStart,newStop);
			filteredMap.add(temp);
		}
		
		return filteredMap;
		
	}
}
