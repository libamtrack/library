package amtrackservice.client;

import amtrackservice.client.Calculation;

/**
 * This class provides the access to the logger windows 
 * to print out warnings and errors.
 * 
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 */
public class Logger {

	private static Calculation c = null;

	/**
	 * sets the MainView on which the logs will be printed.
	 * @param mv
	 */
	public static void init(Calculation c) {
		Logger.c = c;
	}

	/**
	 * adds a new line on the logger windows
	 * @param message the message to be logged
	 */
	public static void log(String message) {
		if( c != null )
			c.addLoggerMessage("[LOG]: " + message);
		System.out.println("[LOG]: " + message);
	}
	
	/**
	 * adds a new line on the logger windows
	 * @param message the message to be logged
	 */
	public static void info(String message) {
		if( c != null )
			c.addLoggerMessage("[INFO]: " + message);
		System.out.println("[INFO]: " + message);
	}
	
	/**
	 * adds a new line on the logger windows
	 * @param message the message to be logged
	 */
	public static void error(String message) {
		if( c != null )
			c.addLoggerMessage("<B>[ERROR]: " + message + "</B>");
		System.out.println("[ERROR]: " + message);
	}
	
	/**
	 * this method will cause logger window to show up
	 */
	public static void show() {
		if( c != null )
			c.showStatus();
	}

	/**
	 * this method will cause logger window to hide
	 */
	public static void hide() {
		if( c != null )
			c.hideStatus();
	}

}
