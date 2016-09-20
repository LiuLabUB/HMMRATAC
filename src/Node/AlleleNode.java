package Node;

public class AlleleNode {
	private String SnpName = "";
	private String ConName = "";
	private int Position = 0;
	private char NT = '\0';
	private int ConStart = 0;
	private int ConStop = 0;
	private String Header = "";
	private StringBuilder Seq;
	private String Quality = "";
	
	public AlleleNode(){
		
	}
	
	public AlleleNode(String snp,String con,int pos,char nt,int start,int stop,String head,StringBuilder seq, String qual){
		SnpName = snp;
		ConName = con;
		Position = pos;
		NT = nt;
		ConStart = start;
		ConStop = stop;
		Header = head;
		Seq = seq;
		Quality = qual;
	}
	public int getPos(){
		return Position;
	}
	public void setPos(int pos){
		Position = pos;
	}
	public int getStart(){
		return ConStart;
	}
	public void setStart(int start){
		ConStart = start;
	}
	public int getStop(){
		return ConStop;
	}
	public void setStop(int stop){
		ConStop = stop;
	}
	
	public String getSnpName(){
		return SnpName;
	}
	public void setSnp(String snp){
		SnpName = snp;
	}
	
	public String getCon(){
		return ConName;
	}
	public void setCon(String con){
		ConName = con;
	}
	
	public char getNT(){
		return NT;
	}
	public void setNT(char nt){
		NT = nt;
	}
	
	public String getHeader(){
		return Header;
	}
	public void setHeader(String head){
		Header = head;
	}
	
	public StringBuilder getSeq(){
		return Seq;
	}
	public void setSeq(StringBuilder seq){
		Seq = seq;
	}
	
	public String getQuality(){
		return Quality;
	}
	public void setQual(String qual){
		Quality = qual;
	}

}
