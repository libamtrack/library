package amtrackservice.client;

import java.util.ArrayList;
import java.util.HashMap;

public class MapList<K,V>{
	private ArrayList<K> keys = null;
	private ArrayList<V> values = null;
	private HashMap<K,V> mapping = null;
	
	public MapList() {
		keys = new ArrayList<K>();
		values = new ArrayList<V>();
		mapping = new HashMap<K,V>();
	}
	
	public void put(K key, V value){
		mapping.put(key, value);
		keys.add(key);
		values.add(value);
	}
	
	public ArrayList<V> getValues(){
		return values;
	}
	
	public ArrayList<K> getKeys(){
		return keys;
	}
	
	public V getValue(K key){
		return mapping.get(key);
	}
}
