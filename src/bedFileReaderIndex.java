package GenomeFileReaders;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;

import Node.TagNode;

public class bedFileReaderIndex {
	
	private File input;
	private File genome;
	
	public bedFileReaderIndex(File i, File g){
		input = i;
		genome = g;
	}
	
	private void read() throws FileNotFoundException{
		GenomeFileReader gReader = new GenomeFileReader(genome);
		
		ArrayList<TagNode> gen = gReader.getMap();
		HashMap<String,ArrayList<TagNode>> IndexHash = new HashMap<String,ArrayList<TagNode>>();
		HashMap<Integer,ArrayList<TagNode>> readHash = new HashMap<Integer,ArrayList<TagNode>>();
		
		int counter = 0;
		for (int i = 0;i < gen.size();i++){
			String chr = gen.get(i).getChrom();
			int stop = gen.get(i).getStop();
			ArrayList<TagNode> tempNode = new ArrayList<TagNode>();
			
			for (int a = 0;a < stop;a += 500000){
				TagNode t = new TagNode(chr,a,a+500000,counter);
				tempNode.add(t);
				counter++;
			}
			IndexHash.put(chr, tempNode);
		}
		
		
		
	}

}
