package WigMath;

import java.io.IOException;
import java.util.ArrayList;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import Node.TagNode;

public class PullWigAboveCutoff {

	private static String input = null;
	private static double cutoff = 0.0;
	private ArrayList<TagNode> g;
	
	private String wig;
	private double cut;
	private ArrayList<TagNode> result;
	public PullWigAboveCutoff(String w,double c, ArrayList<TagNode> G) throws IOException{
		wig = w;
		cut = c;
		g = G;
		set();
	}
	public ArrayList<TagNode> getResults(){return result;}
	private void set() throws IOException{
		result = new ArrayList<TagNode>();
		
		for (int i = 0;i < g.size();i++){
			MeanAndStd results = new MeanAndStd(wig,g.get(i));
			double mean = results.getMean();
			double std = results.getStd();
			results=null;
			BBFileReader reader = new BBFileReader(wig);
			String chrom = g.get(i).getChrom();
			int begin = g.get(i).getStart();
			int end = g.get(i).getStop();
			BigWigIterator iter = reader.getBigWigIterator(chrom,begin,chrom,end,false);
			while (iter.hasNext()){
				WigItem item = iter.next();
				String chr = item.getChromosome();
				int start = item.getStartBase();
				int stop = item.getEndBase();
				double value = item.getWigValue();
				value = (value-mean) / std;
				if (value >= cut){
					TagNode temp = new TagNode(chrom,start,stop,value);
					result.add(temp);
				}
			
			}
		}
		
	}
	
	/*
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			case'i':
				input = args[i+1];
				i++;
				break;
			
			case'c':
				cutoff = Double.parseDouble(args[i+1]);
				i++;
				break;
			case'h':
				printUsage();
				System.exit(1);
			}
		}
		if (input == null || cutoff == 0.0){
			printUsage();
			System.exit(1);
		}
		PullWigAboveCutoff pw = new PullWigAboveCutoff(input,cutoff);
		ArrayList<TagNode> results = pw.getResults();
		for (TagNode r : results){
			System.out.println(r.getChrom()+"\t"+r.getStart()+"\t"+r.getStop()+"\t"+r.getScore2());
		}
		
	}
	*/
	private static void printUsage(){
		System.out.println("Usage: java7 -jar BigWigZScore.jar");
		System.out.println("\tRequired Parameters:");
		System.out.println("\t-i <BigWigFile>\n\t\n\t-c <double> Zscore"
				+ "cutoff to print");
	}
}
