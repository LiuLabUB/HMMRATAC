package Node;
public class FASTQNode {

	
		private String Header = "";
		private StringBuilder Sequence;
		private char Plus = '\0';
		private String Quality = "";
		
		public FASTQNode(){
			
		}
		
		public FASTQNode(String header,StringBuilder seq,char plus,String quality)
		{
			Header = header;
			Sequence = seq;
			Plus = plus;
			Quality = quality;
		}
		
		public String getHeader(){
			return Header;
		}
		public void setHeader(String header){
			Header = header;
		}
		
		public StringBuilder getSeq(){
			return Sequence;
		}
		public void setSeq(StringBuilder seq){
			Sequence = seq;
		}
		
		public char getPlus(){
			return Plus;
		}
		public void setPlus(char plus){
			Plus = plus;
		}
		
		public String getQuality(){
			return Quality;
		}
		public void setQuality(String quality){
			Quality = quality;
		}
}