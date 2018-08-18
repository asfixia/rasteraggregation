
#include <Rcpp.h>
#include <stdlib.h>
#include <algorithm>
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

using namespace Rcpp;



bool isNaN(float value) {
	return value != value;
}

bool isEqual(double a, double b) {
	static double EPSILON = 1E-6;
	return std::abs(a - b) < EPSILON;
}

//Danilo aplica a soma das celulas e retorna o valor. (considerando a % interna da celula).
float focalWeightedSum(double * lines, int xWidth, int yHeight,
	GDALDataType originalMapCellType, double nullValue,
	int fromX, int toX, int fromY, int toY,
	double weightIniLine, double weightToLine, double weightIniCol, double weightToCol) {
	double summed = 0;

	for (int iCurY = fromY; iCurY <= toY; iCurY++) {
		int curIndex = iCurY * xWidth + fromX;

		double weightY;
		if (iCurY == fromY && iCurY == toY)
			weightY = weightIniLine - (1 - weightToLine);
		else if (iCurY == fromY)
			weightY = weightIniLine;
		else if (iCurY == toY)
			weightY = weightToLine;
		else
			weightY = 1;

		for (int iCurX = fromX; iCurX <= toX; iCurX++, curIndex++) {
			//double curValue = getValue(lines, curIndex, originalMapCellType);
			double curValue = lines[curIndex];

			if (!isNaN(curValue) && !isEqual(curValue, nullValue)) {

				double weightX;
				if (iCurX == fromX && iCurX == toX)
					weightX = weightIniCol - (1 - weightToCol);
				else if (iCurX == fromX)
					weightX = weightIniCol;
				else if (iCurX == toX)
					weightX = weightToCol;
				else
					weightX = 1;
				summed += curValue * weightX * weightY;
			}
			//Rcpp::Rcout << "isNan " << isNaN(curValue) << " nulVal: " << nullValue << " isNull: " << isEqual(curValue, nullValue) << " curValue: " << curValue << " summed: " << summed << "\n";
			//Rcpp::Rcout << iCurX << " oriX|oriY " << iCurY << "\n";
		}
	}
	//Rcpp::Rcout << "TotalSummed: " << summed << "\n";

	return summed;
}

float focalWeightedAverage(double * lines, int xWidth, int yHeight,
	GDALDataType originalMapCellType, double nullValue,
	int fromX, int toX, int fromY, int toY,
	double weightIniLine, double weightToLine, double weightIniCol, double weightToCol) {
	double summed = 0;
	double quantity = 0;

	for (int iCurY = fromY; iCurY <= toY; iCurY++) {
		int curIndex = iCurY * xWidth + fromX;

		double weightY;
		if (iCurY == fromY && iCurY == toY)
			weightY = weightIniLine - (1 - weightToLine);
		else if (iCurY == fromY)
			weightY = weightIniLine;
		else if (iCurY == toY)
			weightY = weightToLine;
		else
			weightY = 1;

		for (int iCurX = fromX; iCurX <= toX; iCurX++, curIndex++) {
			double curValue = lines[curIndex];

			if (!isNaN(curValue) && !isEqual(curValue, nullValue)) {
				double weightX;
				if (iCurX == fromX && iCurX == toX)
					weightX = weightIniCol - (1 - weightToCol);
				else if (iCurX == fromX)
					weightX = weightIniCol;
				else if (iCurX == toX)
					weightX = weightToCol;
				else
					weightX = 1;
				quantity += weightX * weightY;
				summed += curValue * weightX * weightY;
			}
			//Rcpp::Rcout << "isNan " << isNaN(curValue) << " nulVal: " << nullValue << " isNull: " << isEqual(curValue, nullValue) << " curValue: " << curValue << " summed: " << summed << "\n";
			//Rcpp::Rcout << iCurX << " oriX|oriY " << iCurY << "\n";

		}
	}
	//Rcpp::Rcout << "TotalSummed: " << summed << "\n";

	return quantity == 0 ? 0 : summed / quantity;
}

enum ResamplingMethod {
	SUM,
	AVERAGE
};


void resampling_algorithm(String originalMapPath, String newMapPath, ResamplingMethod resamplingMethod) {

	GDALDataset  *pOriginalDataset;
	GDALDataset  *pNewDataset;

	//Max Memory settings
	int maxMemoryMB = 1024; //in MB
	int byteToMB = 10 ^ 6;

	GDALAllRegister();

	pOriginalDataset = (GDALDataset *)GDALOpen(originalMapPath.get_cstring(), GA_ReadOnly);
	pNewDataset = (GDALDataset *)GDALOpen(newMapPath.get_cstring(), GA_Update);
	if (pNewDataset == NULL || pOriginalDataset == NULL) {
		exit(1);
	}

	int originalMapXCount = GDALGetRasterXSize(pOriginalDataset);
	int originalMapYCount = GDALGetRasterYSize(pOriginalDataset);

	int newMapXCount = GDALGetRasterXSize(pNewDataset);
	int newMapYCount = GDALGetRasterYSize(pNewDataset);

	GDALRasterBand *originalMapBand = pOriginalDataset->GetRasterBand(1);
	GDALRasterBand *newMapBand = pNewDataset->GetRasterBand(1);

	GDALDataType originalMapCellType = originalMapBand->GetRasterDataType();
	double newMapNullValue = newMapBand->GetNoDataValue();
	double originalMapNullValue = originalMapBand->GetNoDataValue();

	int originalMapFromLineInMemory = INT_MIN;
	int fromLineToWrite = 0;

	//Divido a memória igualmente por enqt.
	int maxSupportedCellByteSize = sizeof(double); // 64 Bytes
	double qntCellsInMemory = maxMemoryMB * byteToMB / maxSupportedCellByteSize;

	double curPropX = double(originalMapXCount) / newMapXCount; //#celulaOsriginal, por celulaNovo
	double curPropY = double(originalMapYCount) / newMapYCount;
	
	int maxMapLinesInMemory = floor(qntCellsInMemory / (originalMapXCount + newMapXCount));
	maxMapLinesInMemory = std::max(maxMapLinesInMemory, (int)ceil(curPropY) + 1);

	double * originalMapLines = (double *)CPLMalloc(sizeof(double) * originalMapXCount * maxMapLinesInMemory);
	double * newMapLines = (double *)CPLMalloc(sizeof(double) * newMapXCount * maxMapLinesInMemory);
	Rcpp::Rcout << curPropX << " X:prop:Y " << curPropY << " " << maxMapLinesInMemory << ":yMem \n";
	//for (int iNewMapY = 0; iNewMapY < 30; iNewMapY++) { //forTests
	for (int iNewMapY = 0; iNewMapY < newMapYCount; iNewMapY++) {
		//Rcpp::Rcout << "\n\niNewY " << iNewMapY << "\n\n\n";
		double realIniLine = curPropY * iNewMapY;
		int iniLine = floor(realIniLine);
		double realToLine = std::min((double)(originalMapYCount - 1), curPropY * (iNewMapY + 1));
		int toLine = floor(realToLine);
		double weightIniLine = 1 - (realIniLine - iniLine);
		double weightToLine = (realToLine - toLine);

		if (toLine >= (originalMapFromLineInMemory + maxMapLinesInMemory - 1)) { //Se não estiver na memória.
			originalMapFromLineInMemory = iniLine;
			int curYSize = std::min(maxMapLinesInMemory, originalMapYCount - iniLine);
			//Rcpp::Rcout << "ReadingOriginalMap: " << iniLine << " F:OriY:T " << toLine << " memOriYSize: " << curYSize << " toOriMemLine: " << (iniLine + curYSize) << " ---\n";
			originalMapBand->RasterIO(GF_Read,
				0, iniLine, originalMapXCount, curYSize, //Map Xoffset, yOffset, xSize, ySize
				originalMapLines, //Data buffer
				originalMapXCount, curYSize, GDT_Float64, //Buffer XSize, ySize, CellType 
				0, 0, 0); //Usarei default.

			/*if (iNewMapY == 28){
				Rcpp::Rcout << "\nRead as array: ";
				for (int iOriY = 0, iCell = 0; iOriY < maxMapLinesInMemory; iOriY++) {
					for (int iOriInd = 0; iOriInd < originalMapXCount; iOriInd++, iCell++)
						Rcpp::Rcout << "(" << iOriInd << ")" << originalMapLines[iCell] << ", ";
					Rcpp::Rcout << "\n";
				}
				Rcpp::Rcout << "\n";
			}*/
		}
		//Rcpp::Rcout << "after read" << "\n";

		for (int iNewMapX = 0; iNewMapX < newMapXCount; iNewMapX++) {
			double realIniCol = curPropX * iNewMapX;
			int iniCol = floor(realIniCol);
			double realToCol = std::min((double)(originalMapXCount - 1), curPropX * (iNewMapX + 1));
			int toCol = floor(realToCol);
			
			double weightIniCol = 1 - (realIniCol - iniCol);
			double weightToCol = (realToCol - toCol); // <- igual // (toCol + 1 - realToCol);
			
			int newMapIndex = (iNewMapY - fromLineToWrite) * newMapXCount + iNewMapX;
			int iniMemLine = iniLine - originalMapFromLineInMemory;
			int toMemLine = toLine - originalMapFromLineInMemory;
			double value;
			switch (resamplingMethod) {
			case SUM:
				value = focalWeightedSum(originalMapLines, originalMapXCount, originalMapYCount,
					GDT_Float64, originalMapNullValue,
					iniCol, toCol, iniMemLine, toMemLine, //Map: fromX, toX, fromMemY, toMemY
					weightIniLine, weightToLine, weightIniCol, weightToCol);
					
				/*if (iNewMapY == 28)
				Rcpp::Rcout << iNewMapX << " X:new:Y " << iNewMapY << " (" << realIniCol << ")" << iniCol << " F:Xori:T " << toCol << "(" << realToCol << ") "
					<< iniLine << " F:Yori:T " << toLine << " " << weightIniCol << " X:WF:Y " << weightIniLine << " " << weightToCol << " X:WT:Y " << weightToLine << " sum:" << value << "---\n";
					*/
				newMapLines[newMapIndex] = isEqual(value, (double)0) ? newMapNullValue : value;
				break;
			case AVERAGE:
				value = focalWeightedAverage(originalMapLines, originalMapXCount, originalMapYCount,
					GDT_Float64, originalMapNullValue,
					iniCol, toCol, iniMemLine, toMemLine, //Map: fromX, toX, fromMemY, toMemY
					weightIniLine, weightToLine, weightIniCol, weightToCol);
				//	Rcpp::Rcout << iNewMapX << " newX:newY " << iNewMapY << " " << iniCol << " oriFX:oriTX " << toCol << " " << iniMemLine << " fOriMY:tOriMY " << toMemLine << " " << iniLine << " oriYF:oriYT " << toLine << " ---\n
				newMapLines[newMapIndex] = isEqual(value, (double)0) ? newMapNullValue : value;
				break;
			default:
				break;
			}
		}

		double curCalculatedYSize = (iNewMapY - fromLineToWrite) + 1;
		//Rcpp::Rcout << iNewMapY << " newMap|max " << maxMapLinesInMemory << " " << curCalculatedYSize << " bufL|" << "\n";
		if (iNewMapY == newMapYCount - 1 || (curCalculatedYSize == maxMapLinesInMemory)) {
			//Rcpp::Rcout << "Writing buff to disk " << fromLineToWrite << " fromLine:toL " << (fromLineToWrite + curCalculatedYSize) << " "  << newMapXCount << " xSize:ySize " << curCalculatedYSize << " ----\n";
			newMapBand->RasterIO(GF_Write,
				0, fromLineToWrite, newMapXCount, curCalculatedYSize, //Map Xoffset, yOffset, xSize, ySize
				newMapLines, //Data buffer
				newMapXCount, curCalculatedYSize, GDT_Float64, //Buffer XSize, ySize, CellType 
				0, 0, 0); //Usarei default.
			fromLineToWrite = iNewMapY + 1;
		}
		//Rcpp::Rcout << "next Y" << "\n";
	}
	CPLFree(originalMapLines);
	CPLFree(newMapLines);
	GDALClose(pOriginalDataset);
	GDALClose(pNewDataset);
}

// [[Rcpp::export]]
void aggregation_resamplingSum(String originalMapPath, String newMapPath) {
	resampling_algorithm(originalMapPath, newMapPath, SUM);
}

// [[Rcpp::export]]
void aggregation_resamplingAverage(String originalMapPath, String newMapPath) {
	resampling_algorithm(originalMapPath, newMapPath, AVERAGE);
}

// [[Rcpp::export]]
String aggregation_gdalVersion() {
	return GDALVersionInfo("VERSION_NUM");
}