package AmTrack;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.security.AccessController;
import java.security.PrivilegedAction;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

public class AmTrackTestGUI extends JPanel 
                             implements ActionListener {

	private static String CONVERT_COMMAND = "Convert";
    private static String EXIT_COMMAND = "Exit";

    private static JTextField fluenceField = new JTextField("1e8", 10);
    private static JTextField sigmaField = new JTextField("2.5", 10);
    private static JTextField nField = new JTextField("0.0", 10);
    private static JTextField fwhmField = new JTextField("0.0", 10);
    
	public native static float[] ATconvertBeamParameters(float fluenceCm2, float sigmaCm2, float N, float fwhmMm);

    public AmTrackTestGUI() {
  //      super(new BorderLayout());
//        panel.setPreferredSize(new Dimension(600, 450));
//        add(panel, BorderLayout.CENTER);
//
//        JPanel panel = new JPanel(new GridLayout(0,2));
//        panel.add(convertButton);
//        panel.add(exitButton);
//	add(panel, BorderLayout.SOUTH);
    	createAndShowGUI(this);
    }

    public void actionPerformed(ActionEvent e) {
        
    	String command = e.getActionCommand();
        if (CONVERT_COMMAND.equals(command)) {
        	String strFluenceCm2 = fluenceField.getText();
        	String strSigmaCm    = sigmaField.getText();
        	String strN          = nField.getText();
        	String strFwhmMm     = fwhmField.getText();
        	float[] results;
        	results = ATconvertBeamParameters(	Float.valueOf(strFluenceCm2).floatValue(), 
							        			Float.valueOf(strSigmaCm).floatValue(), 
							        			Float.valueOf(strN).floatValue(), 
							        			Float.valueOf(strFwhmMm).floatValue());
        	fluenceField.setText(Float.toString(results[0]));
        	sigmaField.setText(Float.toString(results[1]));
        	nField.setText(Float.toString(results[2]));
        	fwhmField.setText(Float.toString(results[3]));
        	System.out.println("Conversion done.");
        } else if (EXIT_COMMAND.equals(command)) {
        	System.out.println("Click the X to exit...");
        }
    }
    

    private static void createAndShowGUI(AmTrackTestGUI gUI) {
    	   
    	// Use priviledged code instead of loadLibrary to call dll from applet
    	
    		try{
    	        AccessController.doPrivileged(new PrivilegedAction(){
    	            public Object run(){
    	            	try
    	            	{
    	            		System.load("D:/workspaces/Eclipse/libamtrack/AT_Java_Release/libAmTrack.dll");
    	            		return null;
    	            	} 
    	                catch (Exception e)
    	                {
    	                	System.out.println("Unable to load libAmTrack.dll");
    	                	return null;
    	                }
    	            }
    	        });
    	   }
           catch (Exception e)
           {
	          	System.out.println("Unable to load libAmTrack.dll");
	       }
    	         
     	//Create and set up the window.
        JFrame frame = new JFrame("libAmTrack Test Applet");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
              
        frame.setLayout(new GridLayout(5,2));
        
        JTextArea fluenceArea = new JTextArea("Fluence / (1/cm2):");
        frame.getContentPane().add(fluenceArea, "1");
        frame.getContentPane().add(fluenceField, "2");
        fluenceField.addActionListener(gUI);
        
        JTextArea sigmaArea = new JTextArea("sigma / cm:");
        frame.getContentPane().add(sigmaArea, "3");
        frame.getContentPane().add(sigmaField, "4");

        JTextArea nArea = new JTextArea("Number of particles:");
        frame.getContentPane().add(nArea, "5");
        frame.getContentPane().add(nField, "6");

        JTextArea fwhmArea = new JTextArea("FWHM / mm:");
        frame.getContentPane().add(fwhmArea, "7");
        frame.getContentPane().add(fwhmField, "8");

        JButton convertButton = new JButton("Convert");
        frame.getContentPane().add(convertButton, "9");
        convertButton.setActionCommand(CONVERT_COMMAND);
        convertButton.addActionListener(gUI);

        JButton exitButton = new JButton("Exit");
        frame.getContentPane().add(exitButton, "10");
        exitButton.setActionCommand(EXIT_COMMAND);
        exitButton.addActionListener(gUI);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
    
    /*
    // called when run as an application
    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }
     */
}
