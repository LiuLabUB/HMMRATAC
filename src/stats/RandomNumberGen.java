package stats;

import java.util.Random;

public class RandomNumberGen {
	
	
	private int _min;
	private int _max;
	
	public RandomNumberGen(int min,int max){
		_min=min;
		_max=max;
	}
	
	public int randInt() {

	    // Usually this can be a field rather than a method variable
	    Random rand = new Random();

	    // nextInt is normally exclusive of the top value,
	    // so add 1 to make it inclusive
	    int randomNum = rand.nextInt((_max - _min) + 1) + _min;

	    return randomNum;
	}

}
