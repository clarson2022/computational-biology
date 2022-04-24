/**
 * Spring 2022
 * Catherine Larson
 * Implementing Expectation-Maximization algorithm for characterizing motifs in genomic sequence.
 */

import java.util.*;
import java.lang.Math;
public class EM_MotifSearch extends MotifSearch {
    
    //constructor
    public EM_MotifSearch(java.lang.String fileName, int motifLength){
        super(fileName, motifLength);
    }

    /**
     * The initial random seed step in the EM algorithm.
     * For each sequence, randomly determine the start index of a motif instance in the sequence.
     */
    public void setRandomLocationsForMotifInstances(){
        for(int i = 0; i<numSequences; i++){
            //get random index
            int randomIndex;
            Random rand = new Random();
            randomIndex = rand.nextInt(sequences.get(i).length()-motifLength+1);
            instanceLocations.set(i, randomIndex);
        }
    }

    /**
     * The Maximization step in the EM algorithm.
     * Based on the motif instances in the sequences, creates a matrix motif model. The resulting matrix model should be updated with pseudocounts so that no entries in the matrix correspond to 0.0.
     */
    public void determineMatrixModel(){
        double[][] newMatrix = new double[motifLength][numSequences];
        //calculate column by column
        int aCount, cCount, gCount, tCount;
        for(int i = 0; i<motifLength;i++){
            //go through each sequence we have
            aCount = cCount = gCount = tCount = 0;
            //have to go through instance locations
            for(int j = 0; j<numSequences; j++){
                //gets char in i'th position
                int motifStart = instanceLocations.get(j);
                char current = sequences.get(j).charAt(motifStart+i);
                if(current == 'A') aCount++;
                else if(current == 'C') cCount++;
                else if(current == 'G') gCount++;
                else if(current == 'T') tCount++;
                
            }
            //calculate percents
            //A's
            matrix[0][i] = (double)aCount / numSequences;
            //C's
            matrix[1][i] = (double)cCount / numSequences;
            //G's
            matrix[2][i]= (double)gCount / numSequences;
            //T's
            matrix[3][i] = (double)tCount / numSequences;
        }

        //add pseudocounts
        addPseudocountsToMatrix();
    }
    
    /**
     * The Expectation step in the EM algorithm. Based on the matrix model, identifies motif 
     * instances in the sequences. One motif instance is identified in each sequence. For each 
     * sequence, the motif instance that best matches the model is chosen.
     */
    public void determineMotifInstances(){
        //for each sequence
        for(int i = 0; i<numSequences; i++){
            //find motif (returns index)
            int location = findMotif(sequences.get(i));
            //put it into the instance locations
            instanceLocations.set(i, location);
        }
    }

    private int findMotif(String s){
        int motifLocation = 0; 
        double bestScore = 0; 
        for(int i = 0; i<s.length()-motifLength+1; i++){
            String current = s.substring(i, i+motifLength);
            //optimize score and choose
            if(getScoreForMotifInstance(current)>bestScore){
                motifLocation = i;
                bestScore = getScoreForMotifInstance(current);
            }
        }
        return motifLocation;
    }

    /**
     * Given a candidate motif instance, returns the score (probability) of that instance based 
     * on the matrix model. The length of the motif instance specified by String s must be 
     * the same as the number of columns in the matrix model.
     */
    public double getScoreForMotifInstance(String s){
        double score = 1;
        for(int i = 0; i < motifLength; i++){
            char current = s.charAt(i);
            if(current == 'A') {
                score = score*matrix[0][i];
            }
            else if(current == 'C'){
                score = score*matrix[1][i];
            }
            else if(current == 'G'){
                score = score*matrix[2][i];
            }
            else if(current == 'T'){
                score = score*matrix[3][i];
            }
        } 
        return score;
    }

    /**
     *  Returns the information content associated with the matrix model.
     */
    public double getInformationContentOfMatrix(){
        //get information content of each position
        double matrixIC = 0;
        for(int i = 0; i<motifLength; i++){
            //add up all IC's of each colum in the matrix
            matrixIC += getInformationContentofPosition(i);
        } 
        return matrixIC;
    }

    /**
     * get information content of matrix at position i
     */
    private double getInformationContentofPosition(int i){
        char[] nucleotides = new char[]{'A', 'C', 'G', 'T'};
        double content = 0;
        double nucContent;
        int count = 0;
        for(char nuc : nucleotides){
            double q = getNucleotideContent(nuc);
            double m = matrix[count][i];
            nucContent = m*(log2(m/q));
            content += nucContent;
            count++;
        }
        return content;
    }

    /**
     * calculates the log base two of a double
     */
    private double log2(double num){
        return Math.log10(num) / Math.log10(2);
    }

    /**
     * Executes the EM (Expectation Maximization) algorithm.
     */
    public void EM(){
        //Initially, the EM algorithm is randomly seeded
        setRandomLocationsForMotifInstances();
        
        boolean converged = false;
        double prevIC = -1;
        
        //Then, the Maximization and Expectation steps are alternately repeated until convergence.
        do{
            determineMatrixModel();
            //determineMotifInstances();
            double currentIC = getInformationContentOfMatrix();
            determineMotifInstances();
            if(currentIC > prevIC){
                //then improving 
                prevIC = currentIC;
            }
            else{
                //not improving
                converged = true;
            }
        }
        while(!converged);
    }

    /**
     * Executes the EM (Expectation Maximization) algorithm multiple times.
     * The number of times that the algorithm is executed is specified by the integer parameter. 
     * The best motif, as determined by information content, is identified over all executions of 
     * the algorithm. Upon completion of this method, this EM_MotifSearch should correspond to the 
     * best motif (including matrix model and motif instances) identified over all executions of the 
     * algorithm.
     */
    public void run_EM_multiple_times(int iterations){
        double bestIC = 0;
        double[][] bestMatrix = new double[4][motifLength];
        Vector<Integer> bestInstanceLocations = new Vector<Integer>();

        for(int i = 0; i<iterations; i++){
            EM();
            if(getInformationContentOfMatrix()>bestIC){
                bestMatrix = getMatrix();
                bestInstanceLocations = getInstanceLocations();
                bestIC = getInformationContentOfMatrix();
            }
        }
        matrix = bestMatrix;
        instanceLocations = bestInstanceLocations;
    }
    

    public static void main(String[] args) {
        EM_MotifSearch em = new EM_MotifSearch(args[0], Integer.parseInt(args[1]));
        em.run_EM_multiple_times(Integer.parseInt(args[2]));

        System.out.println("*****FINAL MOTIFS*****");
        System.out.println(em.matrixToString());
        System.out.println("Motif instances: \n"+ em.motifInstancesToString());
        System.out.println("Consensus Sequence: "+ em.getConsensusSequence());
        System.out.println("Information Content: "+ em.getInformationContentOfMatrix());
    }
}
