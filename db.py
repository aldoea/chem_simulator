import xlrd

sheet = None

def open():
	global sheet
	book=xlrd.open_workbook('db_elements.xlsx')
	sheet=book.sheet_by_index(0)

def getElements():
	global sheet
	components = []
	nrows = sheet.nrows
	for i in range(0,nrows):
		cell = sheet.cell(i,0)
		components.append(cell.value)
	return components	
