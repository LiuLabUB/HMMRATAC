package Node;
public class ORFNode {

	
		private String Chrom = "";
		private int Start = 0;
		private int Stop = 0;
		private int Direction = 0;
		private String ID = "";
		private double FPKM = 0.0;
		
		public ORFNode(){
			
		}
		
		public ORFNode(String chr,int start,int stop,int dir,String name,double fpkm)
		{
			Chrom = chr;
			Start = start;
			Stop = stop;
			Direction = dir;
			ID = name;
			FPKM = fpkm;
		}
		
		public double getFPKM(){
			return FPKM;
		}
		public void setFPKM(double fpkm){
			FPKM = fpkm;
		}
		
		public String getChrom(){
			return Chrom;
		}
		public void setChrom(String chr){
			Chrom = chr;
		}
		
		public int getStart(){
			return Start;
		}
		public void setStart(int start){
			Start = start;
		}
		
		public int getStop(){
			return Stop;
		}
		public void setStop(int stop){
			Stop = stop;
		}
		
		public int getDir(){
			return Direction;
		}
		public void setDir(int dir){
			Direction = dir;
		}
		
		public String getID(){
			return ID;
		}
		public void setID(String name){
			ID = name;
		}

}
