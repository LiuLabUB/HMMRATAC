package WigMath;

import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import Node.TagNode;
import stats.HG19Stats;

public class FoldChange {

	private static String wig = null;
	private static double upper = 10.0;
	private static double lower = 2.0;
	
	private String w = null;
	private double u;
	private double l;
	private ArrayList<TagNode> g;
	private ArrayList<TagNode> result;
	
	public FoldChange(String W,double U,double L,ArrayList<TagNode> G) throws IOException{
		w = W;
		u = U;
		l = L;
		g = G;
		set();
	}
	public ArrayList<TagNode> getResults(){return result;}
	private void set() throws IOException{
		MeanAndStd results = new MeanAndStd(w);
		double mean = results.getMean();
		result = new ArrayList<TagNode>();
		BBFileReader wigReader = new BBFileReader(w);
		for (int i = 0;i < g.size();i++){
			String chrom = g.get(i).getChrom();
			int begin = g.get(i).getStart();
			int end = g.get(i).getStop();
			BigWigIterator iter = wigReader.getBigWigIterator(chrom,begin,chrom,end,false);
			while(iter.hasNext()){
				WigItem item = iter.next();
				String chr = item.getChromosome();
				int start = item.getStartBase();
				int stop = item.getEndBase();
				double value = item.getWigValue();
				if ((value/mean) >= l && (value/mean) <= u){
					TagNode temp = new TagNode(chrom,start,stop,(value/mean));
					result.add(temp);
				}
				
			}
		}
	}
	
	public static void main(String[] args) throws IOException {
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {
			case'w':
				wig = args[i+1];
				i++;
				break;
			case'a':
				upper = Integer.parseInt(args[i+1]);
				i++;
				break;
			case'b':
				lower = Integer.parseInt(args[i+1]);
				i++;
				break;
			case'h':
				printUsage();
				System.exit(1);
			}
		}
		
		if (wig == null){
			printUsage();
			System.exit(1);
		}
		
		
		HG19Stats hg19 = new HG19Stats();
		ArrayList<TagNode> stats = hg19.getStats();
		FoldChange fc = new FoldChange(wig,upper,lower,stats);
		ArrayList<TagNode> results = fc.getResults();
		for (TagNode r : results){
			System.out.println(r.getChrom()+"\t"+r.getStart()+"\t"+r.getStop()+"\t"+r.getScore2());
		}
		
		
		
	}
	private static void printUsage(){
		System.out.println("Usage: java7 -jar BigWigZScore.jar");
		System.out.println("\tRequired Parameters:");
		System.out.println("\t-w <BigWigFile>\n\t\n\t-a <double> Upper Fold Change Limit\n\t-b <double> Lower Fold Change Limit");
	}
}
