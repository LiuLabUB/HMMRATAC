package stats;

import java.awt.Point;

public class RectangularSmooth {
	
	private double[] _smoothed;
	private int _width;
	private int _numPass;
	
	
	public RectangularSmooth(double[] list,int width,int numberOfPasses){
		
		_width = width;
		_numPass = numberOfPasses;
		smooth(list);
	}
	public RectangularSmooth(double[] list, int window){
		
		_numPass = FindPara(window).x;
		_width = FindPara(window).y;
		smooth(list);
	}
	
	public double[] getsmoothedData(){
			return _smoothed;
	}
		
	private void smooth(double[] list){
		
		for (int i = 0;i < _numPass;i++){
				list = onePass(list);
		}
		_smoothed = list;
			
	}
	
	private double[] onePass(double[] list){
		
		double[] tempList = new double[list.length];
		double temp;
		int halfWidth = _width/2;
		for (int i = halfWidth;i < list.length-halfWidth;i++){
			
			temp=list[i];
			for (int a = i+1;a <= i+halfWidth;a++){
				temp+=list[a];
			}
			for (int b = i-halfWidth;b <= i-1;b++){
				temp+=list[b];
			}
			temp = temp/_width;
			

			tempList[i]=temp;
			
		}
		return tempList;
	}
	private  Point FindPara(int window){
		int sqW = (int) Math.sqrt(window);
		int closest = 0;
		Point para = new Point();
		for (int n = sqW;n < window;n++){
			for (int w = sqW;w < window;w++){
				int t = findT(n,w);
				if (isCloser(t,closest,window)){
					closest = t;
					para.x = n;
					para.y = w;
				}
			}
		}
		return para;
	}

	private  int findT(int n,int w){
		return n*w-n+1;
	}
	private  boolean isCloser(int t,int closest,int window){
		int one = Math.abs(window - t);
		int two = Math.abs(window - closest);
		if (one < two){
			return true;
		}
		else{
			return false;
		}
	}
	

}
