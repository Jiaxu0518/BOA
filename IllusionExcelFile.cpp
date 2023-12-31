/******************************************************************************************
Copyright           : 2000-2004, Appache  2.0
FileName            : illusion_excel_file.cpp
Author              : Sail
Version             :
Date Of Creation    : 2009年4月3日
Description         :

Others              :
Function List       :
1.  ......
Modification History:
1.Date  :
Author  :
Modification  :

这个类是从网上下载的，我坐享其成，感谢原来的作者，我只试试是稍稍做了一下修正。
修正包括一些参数的使用不谨慎，bool 改为BOOL等,对于对象关系，我改了一部分，感觉原来的作者对于OO的思路部分不是很清楚。
对于这类东西OLE，我完全不了解，用别人封装的东西感觉还是放心了很多，C++，伟大的C++
http://blog.csdn.net/gyssoft/archive/2007/04/29/1592104.aspx

OLE读写EXCEL都比较慢，所以应该尽量减少OLE的次数
对于读取，还有解决方法，请试用一下预加载的方式，这个方法一次加载所有的读取数据,如此速度就飞快了。
据说写数据是没有什么方法加快的
http://topic.csdn.net/t/20030626/21/1962211.html

增加了一些写入方式的代码，保证可以写入EXCEL数据区，但是对于保存，我发现如果调用CLOSE并且保存的方式，
速度非常慢，我不理解为什么。
所以我吧EXCEL打开了，让你进行后续管理,


******************************************************************************************/




//-----------------------excelfile.cpp----------------

#include "StdAfx.h"
#include "IllusionExcelFile.h"



COleVariant
covTrue((short)TRUE),
covFalse((short)FALSE),
covOptional((long)DISP_E_PARAMNOTFOUND, VT_ERROR);

//
CApplication IllusionExcelFile::excel_application_;


IllusionExcelFile::IllusionExcelFile() :
	already_preload_(FALSE)
{
}

IllusionExcelFile::~IllusionExcelFile()
{
	//
	CloseExcelFile();
}


//初始化EXCEL文件，
BOOL IllusionExcelFile::InitExcel()
{

	//创建Excel 2000服务器(启动Excel) 
	if (!excel_application_.CreateDispatch(_T("Excel.Application"))) //"Excel.Application",NULL
	{
		AfxMessageBox(_T("创建Excel服务失败,你可能没有安装EXCEL，请检查!"));
		return FALSE;
	}

	excel_application_.put_DisplayAlerts(FALSE);
	return TRUE;
}

//
void IllusionExcelFile::ReleaseExcel()
{
	excel_application_.Quit();
	excel_application_.ReleaseDispatch();
	excel_application_ = NULL;
}

//打开excel文件
BOOL IllusionExcelFile::OpenExcelFile(const char *file_name)
{
	CString str(file_name);
	//先关闭
	CloseExcelFile();

	//利用模板文件建立新文档 
	excel_books_.AttachDispatch(excel_application_.get_Workbooks(), true);

	LPDISPATCH lpDis = NULL;
	lpDis = excel_books_.Add(COleVariant(str));
	if (lpDis)
	{
		excel_work_book_.AttachDispatch(lpDis);
		//得到Worksheets 
		excel_sheets_.AttachDispatch(excel_work_book_.get_Worksheets(), true);
		//记录打开的文件名称
		open_excel_file_ = file_name;
		return TRUE;
	}
	else
		return FALSE;
}
BOOL IllusionExcelFile::OpenEmptyExcelFile(const char * file_name)
{
	CString str(file_name);
	//先关闭
	CloseExcelFile();

	//利用模板文件建立新文档 
	excel_books_.AttachDispatch(excel_application_.get_Workbooks(), true);

	LPDISPATCH lpDis = NULL;
	lpDis = excel_books_.Add(vtMissing);
	if (lpDis)
	{
		excel_work_book_.AttachDispatch(lpDis);
		//得到Worksheets 
		excel_sheets_.AttachDispatch(excel_work_book_.get_Worksheets(), true);
		//记录打开的文件名称
		open_excel_file_ = file_name;
		return TRUE;
	}
	else
		return FALSE;
}
//关闭打开的Excel 文件,默认情况不保存文件
void IllusionExcelFile::CloseExcelFile(BOOL if_save)
{
	//如果已经打开，关闭文件
	if (open_excel_file_.IsEmpty() == FALSE)
	{
		//如果保存,交给用户控制,让用户自己存，如果自己SAVE，会出现莫名的等待
		if (if_save)
		{
			ShowInExcel(TRUE);
		}
		else
		{
			//
			excel_work_book_.Close(COleVariant(short(FALSE)), COleVariant(open_excel_file_), covOptional);
			excel_books_.Close();
		}

		//打开文件的名称清空
		open_excel_file_.Empty();
	}



	excel_sheets_.ReleaseDispatch();
	excel_work_sheet_.ReleaseDispatch();
	excel_current_range_.ReleaseDispatch();
	excel_work_book_.ReleaseDispatch();
	excel_books_.ReleaseDispatch();
}

void IllusionExcelFile::SaveasXSLFile(const CString &xls_file)
{
	excel_work_book_.SaveAs(COleVariant(xls_file),
		covOptional,
		covOptional,
		covOptional,
		covOptional,
		covOptional,
		0,
		covOptional,
		covOptional,
		covOptional,
		covOptional,
		covOptional);
	return;
}


int IllusionExcelFile::GetSheetCount()
{
	return excel_sheets_.get_Count();
}


CString IllusionExcelFile::GetSheetName(long table_index)
{
	CWorksheet sheet;
	sheet.AttachDispatch(excel_sheets_.get_Item(COleVariant((long)table_index)), true);
	CString name = sheet.get_Name();
	sheet.ReleaseDispatch();
	return name;
}

//按照序号加载Sheet表格,可以提前加载所有的表格内部数据
BOOL IllusionExcelFile::LoadSheet(long table_index, BOOL pre_load)
{
	LPDISPATCH lpDis = NULL;
	excel_current_range_.ReleaseDispatch();
	excel_work_sheet_.ReleaseDispatch();
	lpDis = excel_sheets_.get_Item(COleVariant((long)table_index));
	if (lpDis)
	{
		excel_work_sheet_.AttachDispatch(lpDis, true);
		excel_current_range_.AttachDispatch(excel_work_sheet_.get_Cells(), true);
	}
	else
	{
		return FALSE;
	}

	already_preload_ = FALSE;
	//如果进行预先加载
	if (pre_load)
	{
		PreLoadSheet();
		already_preload_ = TRUE;
	}

	return TRUE;
}

//按照名称加载Sheet表格,可以提前加载所有的表格内部数据
BOOL IllusionExcelFile::LoadSheet(const char* sheet, BOOL pre_load)
{
	CString str(sheet);
	LPDISPATCH lpDis = NULL;
	excel_current_range_.ReleaseDispatch();
	excel_work_sheet_.ReleaseDispatch();
	lpDis = excel_sheets_.get_Item(COleVariant(str));
	if (lpDis)
	{
		excel_work_sheet_.AttachDispatch(lpDis, true);
		excel_current_range_.AttachDispatch(excel_work_sheet_.get_Cells(), true);

	}
	else
	{
		return FALSE;
	}
	//
	already_preload_ = FALSE;
	//如果进行预先加载
	if (pre_load)
	{
		already_preload_ = TRUE;
		PreLoadSheet();
	}

	return TRUE;
}

//得到列的总数
int IllusionExcelFile::GetColumnCount()
{
	CRange range;
	CRange usedRange;
	usedRange.AttachDispatch(excel_work_sheet_.get_UsedRange(), true);
	range.AttachDispatch(usedRange.get_Columns(), true);
	int count = range.get_Count();
	usedRange.ReleaseDispatch();
	range.ReleaseDispatch();
	return count;
}

//得到行的总数
int IllusionExcelFile::GetRowCount()
{
	CRange range;
	CRange usedRange;
	usedRange.AttachDispatch(excel_work_sheet_.get_UsedRange(), true);
	range.AttachDispatch(usedRange.get_Rows(), true);
	int count = range.get_Count();
	usedRange.ReleaseDispatch();
	range.ReleaseDispatch();
	return count;
}

//检查一个CELL是否是字符串
BOOL IllusionExcelFile::IsCellString(long irow, long icolumn)
{
	CRange range;
	range.AttachDispatch(excel_current_range_.get_Item(COleVariant((long)irow), COleVariant((long)icolumn)).pdispVal, true);
	COleVariant vResult = range.get_Value2();
	//VT_BSTR标示字符串
	if (vResult.vt == VT_BSTR)
	{
		return TRUE;
	}
	return FALSE;
}

//检查一个CELL是否是数值
BOOL IllusionExcelFile::IsCellInt(long irow, long icolumn)
{
	CRange range;
	range.AttachDispatch(excel_current_range_.get_Item(COleVariant((long)irow), COleVariant((long)icolumn)).pdispVal, true);
	COleVariant vResult = range.get_Value2();
	//好像一般都是VT_R8
	if (vResult.vt == VT_INT || vResult.vt == VT_R8)
	{
		return TRUE;
	}
	return FALSE;
}

//
CString IllusionExcelFile::GetCellString(long irow, long icolumn)
{

	COleVariant vResult;
	CString str;
	//字符串
	if (already_preload_ == FALSE)
	{
		CRange range;
		range.AttachDispatch(excel_current_range_.get_Item(COleVariant((long)irow), COleVariant((long)icolumn)).pdispVal, true);
		vResult = range.get_Value2();
		range.ReleaseDispatch();
	}
	//如果数据依据预先加载了
	else
	{
		long read_address[2];
		VARIANT val;
		read_address[0] = irow;
		read_address[1] = icolumn;
		ole_safe_array_.GetElement(read_address, &val);
		vResult = val;
	}

	if (vResult.vt == VT_BSTR)
	{
		str = vResult.bstrVal;
	}
	//整数
	else if (vResult.vt == VT_INT)
	{
		str.Format(_T("%d"), vResult.pintVal);
	}
	//8字节的数字 
	else if (vResult.vt == VT_R8)
	{
		str.Format(_T("%0.0f"), vResult.dblVal);
	}
	//时间格式
	else if (vResult.vt == VT_DATE)
	{
		SYSTEMTIME st;
		VariantTimeToSystemTime(vResult.date, &st);
		CTime tm(st);
		str = tm.Format("%Y-%m-%d");

	}
	//单元格空的
	else if (vResult.vt == VT_EMPTY)
	{
		str = "";
	}

	return str;
}

double IllusionExcelFile::GetCellDouble(long irow, long icolumn)
{
	double rtn_value = 0;
	COleVariant vresult;
	//字符串
	if (already_preload_ == FALSE)
	{
		CRange range;
		range.AttachDispatch(excel_current_range_.get_Item(COleVariant((long)irow), COleVariant((long)icolumn)).pdispVal, true);
		vresult = range.get_Value2();
		range.ReleaseDispatch();
	}
	//如果数据依据预先加载了
	else
	{
		long read_address[2];
		VARIANT val;
		read_address[0] = irow;
		read_address[1] = icolumn;
		ole_safe_array_.GetElement(read_address, &val);
		vresult = val;
	}

	if (vresult.vt == VT_R8)
	{
		rtn_value = vresult.dblVal;
	}

	return rtn_value;
}

//VT_R8
int IllusionExcelFile::GetCellInt(long irow, long icolumn)
{
	int num;
	COleVariant vresult;

	if (already_preload_ == FALSE)
	{
		CRange range;
		range.AttachDispatch(excel_current_range_.get_Item(COleVariant((long)irow), COleVariant((long)icolumn)).pdispVal, true);
		vresult = range.get_Value2();
		range.ReleaseDispatch();
	}
	else
	{
		long read_address[2];
		VARIANT val;
		read_address[0] = irow;
		read_address[1] = icolumn;
		ole_safe_array_.GetElement(read_address, &val);
		vresult = val;
	}
	//
	num = static_cast<int>(vresult.dblVal);

	return num;
}

void IllusionExcelFile::SetCellString(long irow, long icolumn, CString new_string)
{
	COleVariant new_value(new_string);
	CRange start_range = excel_work_sheet_.get_Range(COleVariant(_T("A1")), covOptional);
	CRange write_range = start_range.get_Offset(COleVariant((long)irow - 1), COleVariant((long)icolumn - 1));
	write_range.put_Value2(new_value);
	start_range.ReleaseDispatch();
	write_range.ReleaseDispatch();

}

void IllusionExcelFile::SetCellInt(long irow, long icolumn, int new_int)
{
	COleVariant new_value((long)new_int);

	CRange start_range = excel_work_sheet_.get_Range(COleVariant(_T("A1")), covOptional);
	CRange write_range = start_range.get_Offset(COleVariant((long)irow - 1), COleVariant((long)icolumn - 1));
	write_range.put_Value2(new_value);
	start_range.ReleaseDispatch();
	write_range.ReleaseDispatch();
}

void IllusionExcelFile::SetCellDouble(long irow, long icolumn, double new_double)
{
	COleVariant new_value(new_double);
	CRange start_range = excel_work_sheet_.get_Range(COleVariant(_T("A1")), covOptional);
	CRange write_range = start_range.get_Offset(COleVariant((long)irow - 1), COleVariant((long)icolumn - 1));
	write_range.put_Value2(new_value);
	start_range.ReleaseDispatch();
	write_range.ReleaseDispatch();
}


//
void IllusionExcelFile::ShowInExcel(BOOL bShow)
{
	excel_application_.put_Visible(bShow);
	excel_application_.put_UserControl(bShow);
}

//返回打开的EXCEL文件名称
CString IllusionExcelFile::GetOpenFileName()
{
	return open_excel_file_;
}

//取得打开sheet的名称
CString IllusionExcelFile::GetLoadSheetName()
{
	return excel_work_sheet_.get_Name();
}



//取得列的名称，比如27->AA
char *IllusionExcelFile::GetColumnName(long icolumn)
{
	static char column_name[64];
	size_t str_len = 0;

	while (icolumn > 0)
	{
		int num_data = icolumn % 26;
		icolumn /= 26;
		if (num_data == 0)
		{
			num_data = 26;
			icolumn--;
		}
		column_name[str_len] = (char)((num_data - 1) + 'A');
		str_len++;
	}
	column_name[str_len] = '\0';
	//反转
	_strrev(column_name);

	return column_name;
}

//预先加载
void IllusionExcelFile::PreLoadSheet()
{

	CRange used_range;

	used_range = excel_work_sheet_.get_UsedRange();


	VARIANT ret_ary = used_range.get_Value2();
	if (!(ret_ary.vt & VT_ARRAY))
	{
		return;
	}
	//
	ole_safe_array_.Clear();
	ole_safe_array_.Attach(ret_ary);
}