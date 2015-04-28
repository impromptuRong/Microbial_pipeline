import java.util.*;
import java.io.*;

class select{
	public static void main(String[] args) throws FileNotFoundException {
		Map<String, String[]> fileList = new HashMap<String, String[]>();
		Scanner s = new Scanner(new File("./z.Files.txt"));
		s.useDelimiter("\n");
		while(s.hasNext()){fileList.put(s.next(), new String[2]);}
		
		File folder1 = new File("./AB1ED");
		File folder2 = new File("./AB603");
		File[] LOF = folder2.listFiles();

		for (int i=0; i<LOF.length; i++){
			if(LOF[i].isFile()){
				String filename = LOF[i].getName();
				String[] tmp = filename.split("\\.");
				String type = tmp[tmp.length-1];
				if(type.equals("gz")){
/*					try {
//System.out.println(filename);
						Runtime.getRuntime().exec("gunzip " + "./AB1ED/" + filename);
						Runtime.getRuntime().exec("gunzip " + "./AB603/" + filename);
					} catch (IOException e){
						System.err.println("Error: " + e.getMessage());
					}
*/
					StringBuilder newname = new StringBuilder(tmp[0]);
					for(int j=1; j<tmp.length-1; j++){
						newname.append(".").append(tmp[j]);
					}
					filename = newname.toString();
				}

				tmp = filename.split("\\_");
				String query = tmp[0];
				String region = tmp[tmp.length-2];

				if(fileList.containsKey(query)){
					if(region.equals("R1")){fileList.get(query)[0] = filename;}
					else{fileList.get(query)[1] = filename;}
				}
			}
		}

		try{
			FileWriter fstream = new FileWriter("./stability.files");
			BufferedWriter out = new BufferedWriter(fstream);
			for(String key : fileList.keySet()){
				String f1 = fileList.get(key)[0];
				String f2 = fileList.get(key)[1];
				out.write(key + "\t./fastq/" + f1 + "\t./fastq/" + f2 + "\n");
System.out.println("cat" + " ./AB1ED/" + f1 + " ./AB603/" + f1 + " > ./fastq/" + f1);
System.out.println("cat" + " ./AB1ED/" + f2 + " ./AB603/" + f2 + " > ./fastq/" + f2);
				try {
					Runtime.getRuntime().exec("cat" + " ./AB1ED/" + f1 + " ./AB603/" + f1 + " > ./fastq/" + f1);
					Runtime.getRuntime().exec("cat" + " ./AB1ED/" + f2 + " ./AB603/" + f2 + " > ./fastq/" + f2);
				} catch (IOException e){
					System.err.println("Error: " + e.getMessage());
				}

			}
			out.close();
		} catch (Exception e){
			System.err.println("Error: " + e.getMessage());
		}
		return;
	}
}

