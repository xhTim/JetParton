#include <iostream>
#include <fstream>
#include <vector>
//#include <filesystem>
#include <sstream>
//namespace fs = std::filesystem;

int main() {
    const std::string baseDir = "../../PlaygroundZPC5/";
    const std::string qaqbDir = "Calculate_observables/results/QAQB/";
    const std::string pNchDir = "Calculate_observables/results/pNch/";
    const std::string dNdetaDir = "Calculate_observables/results/dNdeta/";

    double QAQBall[12][82] = {0.0};
    double dNchdeta[91][12] = {0.0};
    double pNch[50][3] = {0.0};

    for (int ifolder = 0; ifolder < 20000; ++ifolder) {
        std::cout << ifolder << std::endl;

        // Check if the QAQB folder exists
        std::string qaqbPath = baseDir + "job-" + std::to_string(ifolder) + "/" + qaqbDir;
        //std::cout << qaqbPath << std::endl;
        /*if (!fs::exists(qaqbPath)) {
            std::cerr << "No such event " << ifolder << std::endl;
            continue;
        }
        */
        std::string fff = qaqbPath  + std::to_string(ifolder);
        // Load QAQBall
        //for (const auto& entry : fs::directory_iterator(qaqbPath)) {
        {
            // Check if the file is open
            std::ifstream inputFile(fff);
            if (!inputFile.is_open()) {
                std::cerr << "Error aopening the file." << std::endl;
                continue; // Return an error code
            }
            
            // Read the data line by line
            std::string line;
            bool firstLine = true;
            int linecount = 0;
            while (std::getline(inputFile, line)) {
               // Skip the first line
               if (firstLine) {
                   firstLine = false;
                   continue;
               }
               // Use a stringstream to parse each line
               std::istringstream iss(line);
               double value;
               int cccount = 0;
               // Read each value from the stringstream
               while (iss >> value) {
                   QAQBall[linecount][cccount] = QAQBall[linecount][cccount] + value;
                   cccount++;
               }
               linecount++;
            }
            // Close the file
            inputFile.close();
        }
    
        // Load pNch
        std::string pNchFile = baseDir + "job-" + std::to_string(ifolder) + "/" + pNchDir;
        /*if (!fs::exists(pNchFile)) {
            std::cerr << "No such event " << ifolder << std::endl;
            continue;
        }
        */
        
        //for (const auto& entry : fs::directory_iterator(pNchFile)) {
        std::string ddd = pNchFile  + std::to_string(ifolder);
        {    // Check if the file is open
            std::ifstream inputFile2(ddd);
            if (!inputFile2.is_open()) {
                std::cerr << "Error bopening the file. " <<  std::endl;
                continue; // Return an error code
            }
            
            // Read the data line by line
            std::string line;
            bool firstLine = true;
            int linecount = 0;
            while (std::getline(inputFile2, line)) {
               // Use a stringstream to parse each line
               std::istringstream isss(line);
               double value;
               int cccount = 0;
               // Read each value from the stringstream
               while (isss >> value) {
                   pNch[linecount][cccount] = pNch[linecount][cccount] + value;
                   cccount++;
               }
               linecount++;
            }
            // Close the file
            inputFile2.close();
        }

        // Load dNchdeta
        std::string dNdetaFile = baseDir + "job-" + std::to_string(ifolder) + "/" + dNdetaDir;
       /* if (!fs::exists(dNdetaFile)) {
            std::cerr << "No such event " << ifolder << std::endl;
            continue;
        }
        */
        //for (const auto& entry : fs::directory_iterator(dNdetaFile)) {
        std::string eee = dNdetaFile  + std::to_string(ifolder);
        {    // Check if the file is open
            std::ifstream inputFile3(eee);
            if (!inputFile3.is_open()) {
                std::cerr << "Error copening the file." << std::endl;
                std::cout << eee << std::endl;
                continue; // Return an error code
            }
            
            // Read the data line by line
            std::string line;
            bool firstLine = true;
            int linecount = 0;
            while (std::getline(inputFile3, line)) {
               // Use a stringstream to parse each line
               std::istringstream issss(line);
               double value;
               int cccount = 0;
               // Read each value from the stringstream
               while (issss >> value) {
                   dNchdeta[linecount][cccount] = dNchdeta[linecount][cccount] + value;
                   cccount++;
               }
               linecount++;
            }
            // Close the file
            inputFile3.close();
        }
    }

    // Save results to files
    std::ofstream qaqbOutput("QAQBall");
    for (int ii=0; ii<12; ii++) {
        for (int jj=0; jj<82; jj++) {
            qaqbOutput << QAQBall[ii][jj] << " ";
        }
        qaqbOutput << std::endl;
    }

    std::ofstream pNchOutput("pNch");
    for (int ii=0; ii<50; ii++) {
        for (int jj=0; jj<3; jj++) {
            pNchOutput << pNch[ii][jj] << " ";
        }
        pNchOutput << std::endl;
    }

    std::ofstream dNchdetaOutput("dNchdeta");
    for (int ii=0; ii<91; ii++) {
        for (int jj=0; jj<12; jj++) {
            dNchdetaOutput << dNchdeta[ii][jj] << " ";
        }
        dNchdetaOutput << std::endl;
    }
    qaqbOutput.close();
    pNchOutput.close();
    dNchdetaOutput.close();
    return 0;
}

