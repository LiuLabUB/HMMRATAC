package Node;

public class StateNode {

	private String chr;
	private int pos;
	private int state;
	
	public StateNode(String c, int p, int s){
		chr = c;
		pos = p;
		state = s;
	}
	public void setChr(String c){
		chr = c;
	}
	public String getChr(){
		return chr;
	}
	public void setPos(int p){
		pos = p;
	}
	public int getPos(){
		return pos;
	}
	public void setState(int s){
		state = s;
	}
	public int getState(){
		return state;
	}
	
}
