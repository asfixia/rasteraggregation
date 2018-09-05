#include <Rcpp.h>
#include <ogrsf_frmts.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <regex>
#include <string>
#include <vector>
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

using namespace Rcpp;

namespace Utils {
	/**
	* Verifica se o valor é NaN.
	* @param value {double} Valor a ser verificado.
	* @return {bool} True se for NaN False caso contrário.
	*/
	bool isNaN(double value) {
		return value != value;
	}

	/**
	* Verifica se dois Doubles são iguais.
	* @param a {double} Verifica se os valores a e b são iguais considerando o erro.
	* @return {bool} Retorna true caso os valores forem iguais, False caso contrário.
	*/
	bool isEqual(double a, double b) {
		static double EPSILON = 1E-6;
		return std::abs(a - b) < EPSILON;
	}
	
	/**
	 * Divide a srting str, pelo caractere split, remove os caracteres split.
	 * @param {std::string} String a ser modificada.
	 * @param {char} Caractere que causará split.
	 * @return {std::vector<std::string>} Array de strings separadas pelo split.
	 */
	std::vector<std::string> splitString(std::string str, char split) {
		std::vector<std::string> strings;
		std::istringstream f(str);
		std::string auxStr;
		while (getline(f, auxStr, split)) {
			strings.push_back(auxStr);
		}
		return strings;
	}

	/**
	 * Define um range dos caracteres aceitos.
	 * @param {char} Caractere a ser avaliado.
	 * @return {bool} True caso o caractere seja válido, False caso contrário.
	 */
	bool invalidChar(char c) {
		return !(c >= 0 && c <128);
	}
	/**
	 * Remove os caracteres unicode da string.
	 * @param {std::string &} String que terá os caracteres não unicode removidos.
	 */ 
	void stripUnicode(std::string & str) {
		str.erase(remove_if(str.begin(), str.end(), invalidChar), str.end());
	}

	/**
	 * Faz parse da string do filter, um CSV contendo os valores e operações na forma polonesa inversa ('VALOR',OPERAÇÃO).
	 * @param {std::string} Conteúdo do filtro, em formato CSV, contendo valores e operações. (interpretado por SLD)
	 * @return {std::vector<std::string>} Retorna um array das operações e valores.
	 */
	std::vector<std::string> getParamsFromFilter(std::string curFilter) {
		std::string paramStr = "";
		stripUnicode(curFilter);
		std::vector<std::string> parsedLine = Utils::splitString(curFilter, ',');
		std::vector<std::string> params;
		for (int iStr = 0, strCount = parsedLine.size(); iStr < strCount; iStr++) {
			paramStr += parsedLine[iStr];
			if (paramStr.length() < 2 || paramStr.at(paramStr.length() - 2) != '\\' || iStr + 1 == strCount) {
				//std::cout << paramStr << std::endl;
				params.push_back(paramStr);
				paramStr = "";
			} else {
				paramStr += ",";
			}
		}
		return params;
	}

	
	/**
	 * Tipos de dados que pode-se armazenar no ValueExp.
	 */
	enum ValueType {
		TSTRING,
		TBOOL,
		TDOUBLE
	};

	class ValueExp {
	public:
		std::string strValue = "";
		bool boolValue = false;
		double doubleValue = 0;
		ValueType type;

		ValueExp(std::string value) {
			strValue = value;
			type = ValueType::TSTRING;
		}

		ValueExp(bool value) {
			boolValue = value;
			type = ValueType::TBOOL;
		}

		ValueExp(double value) {
			doubleValue = value;
			type = ValueType::TDOUBLE;
		}

		bool andValue(ValueExp v) {
			if (type == v.type && type == TBOOL)
				return boolValue && v.boolValue;
			else
				return false;
		}

		int compare(ValueExp v) {
			int result;
			if (v.type != type)
				result = false;
			else
				switch (type) {
				case TSTRING:
					result = strValue.compare(v.strValue);
					break;
				case TBOOL:
					result = boolValue - v.boolValue;
					break;
				case TDOUBLE:
					result = std::abs(doubleValue - v.doubleValue) < 0.00001 ? 0 : doubleValue - v.doubleValue;
					break;
				default:
					result = false;
				}
			return result;
		}

		bool orValue(ValueExp v) {
			if (type == v.type && type == TBOOL)
				return boolValue || v.boolValue;
			else
				return false;
		}

		/**
		* Compara sem usar regex, mas aceita * e . para aceitar qualquer caractere.
		*/
		bool like(ValueExp v) {
			if (type == v.type && type == TSTRING) {
				bool isEqual = strValue.length() == v.strValue.length();
				for (int iChar = 0, charCount = strValue.length(); isEqual && iChar < charCount; iChar++) {
					isEqual = strValue.at(iChar) == v.strValue.at(iChar) || strValue.at(iChar) == '*' || strValue.at(iChar) == '.';
				}
				return isEqual; //std::regex_match(strCopy, v.strValue);
			}
			else {
				return false;
			}
		}

		bool result() {
			int result;
			switch (type) {
			case TSTRING:
				result = strValue.empty() == false;
				break;
			case TBOOL:
				result = boolValue;
				break;
			case TDOUBLE:
				result = doubleValue;
				break;
			default:
				result = false;
			}
			return result;
		}
	};

	/**
	 * Pega o objeto no topo o removendo.
	 * @param {std::vector<ValueExp&>} Vetor onde os elementos estão armazenados.
	 * @return {ValueExp} Retorna o elemento do topo.
	 */
	ValueExp getRemoveTop(std::vector<ValueExp>& vector) {
		ValueExp curValue = vector.back();
		vector.pop_back();
		return curValue;
	}
}


/**
* Metodo de resampling a ser aplicado.
*/
enum ResamplingMethod {
	SUM,
	AVERAGE
};

/**
* Soma os elementos da matriz na janela. (Considerando a % interna de sobreposição da célula)
*
* @param {double *} lines Array de valores, representando os dados matriciais.
* @param {int} xWidth Quantidade de colunas da matriz.
* @param {int} yHeight Quantidade de linhas da matriz.
* @param {double} nullValue Valor do nulo ou "sem Informação".
* @param {int} fromX Define a coluna inicial da janela de soma, indexada em 0. (se der from 0, inclui a linha 0)
* @param {int} fromY Define a linha inicial da janela de soma, indexada em 0 (se der to 0, inclui a linha 0)
* @param {int} toX Define a coluna final da janela de soma, indexada em 0. (se der from 0, inclui a linha 0)
* @param {int} toY Define a linha final da janela de soma, indexada em 0 (se der to 0, inclui a linha 0)
* @param {double} weightIniLine Peso atribuido a primeira linha da janela.
* @param {double} weightToLine Peso atribuido a ultima linha da janela.
* @param {double} weightIniCol Peso atribuido a primeira coluna da janela.
* @param {double} weightToCol Peso atribuido a ultima coluna da janela.
* @return {double} Valor da soma dos elementos da janela com os pesos dados.
*/
double focalWeightedSum(double * lines, int xWidth, int yHeight,
	double nullValue, int fromX, int toX, int fromY, int toY,
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

			if (!Utils::isNaN(curValue) && !Utils::isEqual(curValue, nullValue)) {

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

/**
* Calcula a média dos elementos da janela.
*
* @param {double *} lines Array de valores, representando os dados matriciais.
* @param {int} xWidth Quantidade de colunas da matriz.
* @param {int} yHeight Quantidade de linhas da matriz.
* @param {double} nullValue Valor do nulo ou "sem Informação".
* @param {int} fromX Define a coluna inicial da janela de soma, indexada em 0. (se der from 0, inclui a linha 0)
* @param {int} fromY Define a linha inicial da janela de soma, indexada em 0 (se der to 0, inclui a linha 0)
* @param {int} toX Define a coluna final da janela de soma, indexada em 0. (se der from 0, inclui a linha 0)
* @param {int} toY Define a linha final da janela de soma, indexada em 0 (se der to 0, inclui a linha 0)
* @param {double} weightIniLine Peso atribuido a primeira linha da janela.
* @param {double} weightToLine Peso atribuido a ultima linha da janela.
* @param {double} weightIniCol Peso atribuido a primeira coluna da janela.
* @param {double} weightToCol Peso atribuido a ultima coluna da janela.
* @return {double} Valor da média dos elementos na janela
*/
double focalWeightedAverage(double * lines, int xWidth, int yHeight,
	double nullValue, int fromX, int toX, int fromY, int toY,
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

			if (!Utils::isNaN(curValue) && !Utils::isEqual(curValue, nullValue)) {
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


/**
* Calcula a janela necessária e aplica a operação respectiva ao metodo desejado.
* Essa função é responsável por ler e escrever os dados.
*
* PS: O newMapPath é considerado como o mapa já criado, com os tamanhos de célula corretos e ancorado na mesma posição do mapa originalMapPath.
*
* @param originalMapPath {String} Caminho do mapa original.
* @param newMapPath {String} Caminho do mapa a ser reamostrado.
* @param resamplingMethod {ResamplingMethod} Metodo de sampling a ser aplicado.
*/
void resampling_algorithm(String originalMapPath, String newMapPath, ResamplingMethod resamplingMethod) {

	GDALDataset  *pOriginalDataset;
	GDALDataset  *pNewDataset;

	//Max Memory settings
	int maxMemoryMB = 1024; //in MB
	int byteToMB = 10 ^ 6;

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

	//GDALDataType originalMapCellType = originalMapBand->GetRasterDataType();
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
					originalMapNullValue, iniCol, toCol, iniMemLine, toMemLine, //Map: nullValue(noData), fromX, toX, fromMemY, toMemY
					weightIniLine, weightToLine, weightIniCol, weightToCol);

				/*if (iNewMapY == 28)
				Rcpp::Rcout << iNewMapX << " X:new:Y " << iNewMapY << " (" << realIniCol << ")" << iniCol << " F:Xori:T " << toCol << "(" << realToCol << ") "
				<< iniLine << " F:Yori:T " << toLine << " " << weightIniCol << " X:WF:Y " << weightIniLine << " " << weightToCol << " X:WT:Y " << weightToLine << " sum:" << value << "---\n";
				*/
				newMapLines[newMapIndex] = Utils::isEqual(value, (double)0) ? newMapNullValue : value;
				break;
			case AVERAGE:
				value = focalWeightedAverage(originalMapLines, originalMapXCount, originalMapYCount,
					originalMapNullValue, iniCol, toCol, iniMemLine, toMemLine, //Map: nullValue(noData), fromX, toX, fromMemY, toMemY
					weightIniLine, weightToLine, weightIniCol, weightToCol);
				//	Rcpp::Rcout << iNewMapX << " newX:newY " << iNewMapY << " " << iniCol << " oriFX:oriTX " << toCol << " " << iniMemLine << " fOriMY:tOriMY " << toMemLine << " " << iniLine << " oriYF:oriYT " << toLine << " ---\n
				newMapLines[newMapIndex] = Utils::isEqual(value, (double)0) ? newMapNullValue : value;
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

/**
* Aplica o metodo de resampling por soma ponderada do mapa original armazenando o resultado no newMap.
*
* PS: O newMap precisa existir e ter a mesma projeção do mapa originalMapPath.
*
* @param originalMapPath {String} Caminho do mapa original.
* @param newMapPath {String} Caminho do mapa a ser reamostrado.
*/
// [[Rcpp::export]]
void aggregation_resamplingSum(String originalMapPath, String newMapPath) {
	resampling_algorithm(originalMapPath, newMapPath, SUM);
}

/**
* Aplica o metodo de resampling por média ponderada do mapa original armazenando o resultado no newMap.
*
* PS: O newMap precisa existir e ter a mesma projeção do mapa originalMapPath.
*
* @param originalMapPath {String} Caminho do mapa original.
* @param newMapPath {String} Caminho do mapa a ser reamostrado.
*/
// [[Rcpp::export]]
void aggregation_resamplingAverage(String originalMapPath, String newMapPath) {
	resampling_algorithm(originalMapPath, newMapPath, AVERAGE);
}

/**
 * Pega todos os valores unicos do shapefile, no campo 'fieldName'.
 * @param shapePath {Rcpp::String} Caminho do arquivo de shape.
 * @param fieldName {Rcpp::String} Nome do campo do shape.
 * @return {CharacterVector} Array de strings, que cada um representa um valor unico do mapa.
 */
// [[Rcpp::export]]
CharacterVector aggregation_fieldUniqueValuesString(String shapePath, String fieldName) {
	GDALDataset       *poDS =
		(GDALDataset*)GDALOpenEx(shapePath.get_cstring(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (poDS == NULL) {
		printf("Open failed.\n");
		exit(1);
	}

	std::set<std::string> uniqueValues;
	OGRLayer * resultLayer = poDS->ExecuteSQL((
		std::string("SELECT \"") + std::string(fieldName.get_cstring())
		+ std::string("\"  FROM geology")
		).c_str(),
		nullptr,
		nullptr);
	for (OGRFeature * poFeature = resultLayer->GetNextFeature(); poFeature != nullptr;
		poFeature = resultLayer->GetNextFeature()) {
		uniqueValues.insert(poFeature->GetFieldAsString(0));
		poFeature->DestroyFeature(poFeature);
	}
	poDS->ReleaseResultSet(resultLayer);

	int n = uniqueValues.size();
	CharacterVector out(n);

	int i = 0;
	for (auto it = uniqueValues.begin(); it != uniqueValues.end(); ++it, ++i) {
		out[i] = *it;
	}

	return out;
}

/**
 * Soma a área das propriedades que respeitam uma clausula 'where'.
 * @param shapePath {Rcpp::String} Caminho do arquivo shapefile.
 * @param where {Rcpp::String} Clausula where do OGRSql.
 * @return {Double} Soma a área dos atributos.
 */
// [[Rcpp::export]]
double aggregation_sumAreaWhere(String shapePath, String where) {
	GDALDataset       *poDS =
		(GDALDataset*)GDALOpenEx(shapePath.get_cstring(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (poDS == NULL) {
		printf("Open failed.\n");
		exit(1);
	}
	OGRLayer  * resultLayer = poDS->ExecuteSQL(
		(
			std::string("SELECT SUM(OGR_GEOM_AREA) FROM geology WHERE ") + std::string(where.get_cstring())
			).c_str(),
		nullptr,
		nullptr);

	double totalArea = 0;
	for (OGRFeature * poFeature = resultLayer->GetNextFeature(); poFeature != nullptr;
		poFeature = resultLayer->GetNextFeature()) {
		totalArea += poFeature->GetFieldAsDouble(0);
		poFeature->DestroyFeature(poFeature);
	}
	poDS->ReleaseResultSet(resultLayer);
	return totalArea;
}

/**
 * Soma a área dos atributos que tem o 'fieldName' igual á 'fieldValue'.
 * @param shapePath {Rcpp::String} Caminho do shapefile.
 * @param fieldName {Rcpp::String} Nome do campo a ser interpretado.
 * @param fieldValue {Rcpp::String} Valor igual a ser filtrado.
 * @return Retorna a soma dos atributos que tem o 'fieldValue'.
 */
// [[Rcpp::export]]
double aggregation_sumAreaByValue(String shapePath, String fieldName, String fieldValue) {
	String where = std::string("\"") + std::string(fieldName.get_cstring()) + std::string("\"= '")
		+ std::string(fieldValue.get_cstring()) + std::string("'");
	return aggregation_sumAreaWhere(shapePath, where);
}

/**
* Função que retorna uma tabela da area da região agrupada pela coluna.
* Para cada entrada única, são somadas a áreas das features que tem aquele valor.
*
* @param {String} shapePath Caminho do Shapefile.
* @param {String} fieldName Nome da coluna para agrupar.
* @return {Rcpp::DataFrame} Retorna as tuplas [<Propriedade> -> <Area>],..
*/
// [[Rcpp::export]]
DataFrame aggregation_sumAreaGroupedByColumn(String shapePath, String fieldName) {
	CharacterVector allValues = aggregation_fieldUniqueValuesString(shapePath, fieldName);
	NumericVector allAreas(allValues.length());

	auto itArea = allAreas.begin();
	for (CharacterVector::iterator itChar = allValues.begin(); itChar != allValues.end(); ++itChar, ++itArea) {
		*itArea = aggregation_sumAreaByValue(shapePath, fieldName, *itChar);
	}
	return DataFrame::create(Named("Name") = allValues, Named("Area") = allAreas, Named("stringsAsFactors") = false);
}

/**
* Função que retorna a versão do gdal.
* @return {String} Versão do gdal.
*/
// [[Rcpp::export]]
String aggregation_init() {
	GDALAllRegister();
	return GDALVersionInfo("VERSION_NUM");
}

bool interporetExpression(std::string curFilter, Utils::ValueExp ruleValueExp) {
	//Rcpp::Rcout << std::string("curFilter: ") << std::string(curFilter) << std::string("\n");
	Utils::stripUnicode(curFilter);
	if (ruleValueExp.type == Utils::ValueType::TSTRING) {
		//std::cout << std::string("ruleValue: ") << std::string(ruleValueExp.strValue) << std::string("\n");
		Utils::stripUnicode(ruleValueExp.strValue);
	}
	std::vector<Utils::ValueExp> curStack;
	for (std::string curString : Utils::getParamsFromFilter(curFilter)) {
		//std::cout << "curInterpret: " << curString << std::endl;
		if (curString.at(0) == '\'' && curString.at(curString.length() - 1) == '\'') {
			if (ruleValueExp.type == Utils::ValueType::TDOUBLE)
				curStack.push_back(Utils::ValueExp(stod(curString.substr(1, curString.length() - 2))));
			else
				curStack.push_back(Utils::ValueExp(curString.substr(1, curString.length() - 2)));
		} else if (curString.compare("==") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp otherValue = ruleValueExp;
			curStack.push_back(otherValue.compare(value) == 0);
		} else if (curString.compare("!=") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp otherValue = ruleValueExp;
			curStack.push_back(otherValue.compare(value) != 0);
		} else if (curString.compare("<") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp otherValue = ruleValueExp;
			curStack.push_back(otherValue.compare(value) < 0);
		} else if (curString.compare("<=") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp otherValue = ruleValueExp;
			curStack.push_back(otherValue.compare(value) <= 0);
		} else if (curString.compare(">") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp otherValue = ruleValueExp;
			curStack.push_back(otherValue.compare(value) > 0);
		} else if (curString.compare(">=") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp otherValue = ruleValueExp;
			curStack.push_back(otherValue.compare(value) >= 0);
		} else if (curString.compare("LIKE") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp otherValue = ruleValueExp;
			curStack.push_back(otherValue.like(value));
		} else if (curString.compare("NOT") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			curStack.push_back(Utils::ValueExp(value.boolValue));
		} else if (curString.compare("&&") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp other = Utils::getRemoveTop(curStack);
			curStack.push_back(Utils::ValueExp(other.andValue(value)));
		} else if (curString.compare("||") == 0) {
			Utils::ValueExp value = Utils::getRemoveTop(curStack);
			Utils::ValueExp other = Utils::getRemoveTop(curStack);
			curStack.push_back(Utils::ValueExp(other.orValue(value)));
		}
	}
	return (curStack.empty() ? ruleValueExp : Utils::getRemoveTop(curStack)).result();
}

// [[Rcpp::export]]
bool aggregation_interpretExpressionDouble(String curFilter, double ruleValue) {
	return interporetExpression(std::string(curFilter.get_cstring()), Utils::ValueExp(ruleValue));
}

// [[Rcpp::export]]
bool aggregation_interpretExpressionString(String curFilter, String ruleValue) {
	return interporetExpression(std::string(curFilter.get_cstring()), Utils::ValueExp(std::string(ruleValue.get_cstring())));
}
