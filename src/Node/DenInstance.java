package Node;

import net.sf.javaml.core.DenseInstance;

public class DenInstance extends DenseInstance{

	private String _id;
	public DenInstance(double[] att,String ID) {
		super(att);
		_id = ID;
	}
	
	
	public String getStringID(){
	
		return _id;
	}

}
