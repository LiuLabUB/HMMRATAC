package Node;

public class TemplateNode {
	
	private int _chrom;
	private int _pos;
	private char _dir;
	private int _count;
	
	public TemplateNode(int c, int p, char d, int ct){
		_chrom = c;
		_pos = p;
		_dir = d;
		_count = ct;
	}

	public void setChrom(int chr){_chrom = chr;}
	public int getChrom(){return _chrom;}
	public void setPos(int p){_pos = p;}
	public int getPos(){return _pos;}
	public void setDir(char d){_dir = d;}
	public char getDir(){return _dir;}
	public void setCount(int count){_count = count;}
	public int getCount(){return _count;}
	
}
