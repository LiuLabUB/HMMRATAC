package Node;

public class OverlapNode {

	
	private TagNode _node;
	private boolean _hit;
	private boolean _isConsumed;
	
	public OverlapNode(TagNode n,boolean h,boolean c){
		_node = n;
		_hit = h;
		_isConsumed = c;
	}
	public boolean isConsumed(){return _isConsumed;}
	public void setConsumed(boolean c){_isConsumed = c;}
	public boolean hasHit(){return _hit;}
	public TagNode getHit(){return _node;}
	public void setHit(boolean h){_hit = h;}
	public void setOverlap(TagNode n){_node = n;}
}
