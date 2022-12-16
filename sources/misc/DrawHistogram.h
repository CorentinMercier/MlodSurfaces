// --------------------------------------------------------------------------
// This file is part of the reference implementation for the paper
//    Moving Level-of-Detail Surfaces.
//    C. Mercier, T. Lescoat, P. Roussillon, T. Boubekeur, and J-M. Thiery
//    ACM Transaction On Graphics 2022
//    DOI: 10.1145/3528223.3530151
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.
// --------------------------------------------------------------------------

#ifndef DRAWHISTOGRAM_H
#define DRAWHISTOGRAM_H

#include <QGraphicsScene>

class DrawRepartition{
	QGraphicsScene * scene;

public:
	DrawRepartition()
	{
		scene = new QGraphicsScene;
	}
	void clear()
	{
		scene->clear();
		scene->setBackgroundBrush( QBrush(QColor(255,255,255)) );
	}

	void buildSceneFromData( std::vector< double > & repartition , double minvalue , double maxvalue , double w , double h )
	{
		for( unsigned int r = 0 ; r < repartition.size() ; ++r )
		{
			float x_min = (float)(r) / (float)(repartition.size());
			float x_max = (float)(r + 1) / (float)(repartition.size());
			float y_max = ((float)(repartition[r]) - minvalue) / (float)(maxvalue - minvalue);

			QPen pen(QColor(   (int)((1-y_max) * 255) , 0 ,  (int)((y_max) * 255) ));
			QBrush brush(QColor(   (int)((1-y_max) * 255) , 0 ,  (int)((y_max) * 255) ));
			brush.setStyle(Qt::SolidPattern);

			scene->addRect( QRectF(w*x_min, h*(1-y_max) , w*(x_max-x_min) , h*y_max)  ,  pen , brush);
		}
	}

	void buildSceneFromData( std::vector< double > & repartition , double minvalue  , double w , double h )
	{
		double maxvalue;
		for( unsigned int r = 0 ; r < repartition.size() ; ++r )
			maxvalue = std::max( maxvalue , repartition[r] );

		buildSceneFromData(repartition,minvalue,maxvalue,w,h);
	}

	void buildSceneFromData( const std::string & filename , double minvalue  , double w , double h )
	{
		std::vector< double > repartition;

		std::ifstream myfile;
		myfile.open(filename.c_str());
		if (!myfile.is_open())
		{
			std::cout << filename << " cannot be opened" << std::endl;
			return;
		}

		int nData;

		myfile >> nData;

		for(int i = 0 ; i < nData ; ++i )
		{
			double data;
			myfile >> data;
			repartition.push_back(data);
		}

		myfile.close();

		buildSceneFromData(repartition,minvalue,w,h);
	}

	void compareRescaledData(  const std::string & filename1 , const std::string & filename2 , double minvalue  , double w , double h   )
	{
		std::vector< double > repartition1 , repartition2;
		int nData;
		std::ifstream myfile;

		myfile.open(filename1.c_str());
		if (!myfile.is_open())
		{
			std::cout << filename1 << " cannot be opened" << std::endl;
			return;
		}
		myfile >> nData;
		for(int i = 0 ; i < nData ; ++i )
		{
			double data;
			myfile >> data;
			repartition1.push_back(data);
		}
		myfile.close();

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		myfile.open(filename2.c_str());
		if (!myfile.is_open())
		{
			std::cout << filename2 << " cannot be opened" << std::endl;
			return;
		}
		myfile >> nData;
		assert( nData == (int)( repartition1.size() )   &&   "You can only compare tables with same length" );
		for(int i = 0 ; i < nData ; ++i )
		{
			double data;
			myfile >> data;
			repartition2.push_back(data);
		}
		myfile.close();

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		// in our case:
		double max1 = repartition1[0] , max2 = repartition2[0];

		// we create repartition2 / repartition1:

		std::vector< double > repartition;
		buildSceneFromData(repartition,minvalue,w,h);
		for( unsigned int i = 0 ; i < repartition1.size() ; ++i )
			repartition.push_back( max1 * repartition2[i] / ( max2 * repartition1[i] ) );

		buildSceneFromData(repartition,minvalue,w,h);
	}

	void renderImage( const std::string & filename , int w , int h )
	{
		QImage _image(w,h  ,  QImage::Format_RGB32);

		QPainter painter(&_image);
		scene->render(&painter);

		_image.save(QString::fromStdString(filename));
	}
};




class DrawHistogram{
	// use a QGraphicScene to render in the QPainter of a QImage, then save it

	/* // From http://doc.qt.nokia.com/latest/qgraphicsscene.html:
QGraphicsScene scene;
 scene.addItem(...
 ...
 QPrinter printer(QPrinter::HighResolution);
 printer.setPaperSize(QPrinter::A4);

 QPainter painter(&printer);
 scene.render(&painter);
*/
	QGraphicsScene * scene;
	std::vector< int > repartition;
	int repartition_max;
	int n_data;

public:
	DrawHistogram()
	{
		scene = new QGraphicsScene;
	}
	void clear()
	{
		scene->clear();
		scene->setBackgroundBrush( QBrush(QColor(255,255,255)) );
	}


	template< typename T >
	void fillData( const std::vector< T > & data , T data_min_value , T data_max_value , int n_intervalles )
	{
		repartition.clear();
		n_data = 0;
		repartition_max = -1;
		repartition.resize( n_intervalles , 0 );
		for( unsigned int d = 0 ; d < data.size() ; ++d )
		{
			if(  data[d]  >  data_max_value  )
				continue;
			if(  data[d]  <  data_min_value  )
				continue;
			int i_intervalle = std::min( n_intervalles - 1 , (int)( n_intervalles * (data[d] - data_min_value)/(data_max_value - data_min_value) ) );
			++repartition[i_intervalle];
			repartition_max = std::max(  repartition_max  ,  repartition[i_intervalle]  );
			++n_data;
		}
	}

	template< typename T >
	void fillData( const std::vector< T > & data , int n_intervalles = 512 )
	{
		assert( data.size() > 0 );
		T data_max_value = data[0];
		T data_min_value = data[0];
		for( T d : data ) {
			data_max_value = std::max< T > ( data_max_value , d );
			data_min_value = std::min< T > ( data_min_value , d );
		}
		std::cout << "Min in data : " << data_min_value << "  |  Max in data : " << data_max_value << "   |  (" << data.size() << " data points)" << std::endl;
		fillData( data , data_min_value , data_max_value , n_intervalles );
	}


	template< typename T >
	void fillData( const std::vector< T > & data , const std::vector< T > & dataSteps , T data_min_value , T data_max_value )
	{
		unsigned int n_intervalles = dataSteps.size() - 1;
		repartition.clear();
		n_data = 0;
		repartition_max = -1;
		repartition.resize( n_intervalles , 0 );
		for( unsigned int d = 0 ; d < data.size() ; ++d )
		{
			if(  data[d]  >  data_max_value  )
				continue;
			if(  data[d]  <  data_min_value  )
				continue;
			unsigned int i_intervalle = 0;
			for( unsigned int int_it = 0 ; int_it < n_intervalles ; ++int_it )
			{
				i_intervalle = int_it;
				if( data[d] > dataSteps[int_it] )
					break;
			}
			++repartition[i_intervalle];
			repartition_max = std::max(  repartition_max  ,  repartition[i_intervalle]  );
			++n_data;
		}
	}


	template< typename T >
	void fillData(
			const std::vector< T > & data , T data_min_value , T data_max_value , float relevant , int n_intervalles,
			bool keepMinValue = false , bool keepMaxValue = false )
	{
		std::vector< T > values;
		for( unsigned int d = 0 ; d < data.size() ; ++d )
		{
			if(  data[d]  >  data_max_value  )
				continue;
			if(  data[d]  <  data_min_value  )
				continue;
			values.push_back(data[d]);
		}
		std::sort(values.begin(),values.end() , std::less<T>());
		int n_relevant = std::max( 0 , std::min(  (int)( (float)(values.size()) * relevant ) ,  (int)values.size() ) );

		// find the n_relevant values that gives the least range:
		int d_relevant = 0;
		T least_range = std::abs( values[ 0 + n_relevant - 1 ] - values[ 0 ] );
		for( int d = 0 ; d < (int)data.size() - n_relevant ; ++d )
		{
			T range = std::abs( values[ d + n_relevant - 1 ] - values[ d ] );
			if( range < least_range )
			{
				least_range = range;
				d_relevant = d;
			}
		}
		if( !keepMinValue )
			data_min_value = values[d_relevant];
		if( !keepMaxValue )
			data_max_value = values[d_relevant + n_relevant - 1];

		std::cout << "found " << n_relevant << " relevant values between " << data_min_value << " and " << data_max_value << std::endl;

		repartition.clear();
		n_data = 0;
		repartition_max = -1;
		repartition.resize( n_intervalles , 0 );


		for( int d = d_relevant ; d < d_relevant + n_relevant ; ++d )
		{
			int i_intervalle = std::max( 0 , std::min( n_intervalles - 1 , (int)( n_intervalles * (values[d] - data_min_value)/(data_max_value - data_min_value) ) ) );
			++repartition[i_intervalle];
			repartition_max = std::max(  repartition_max  ,  repartition[i_intervalle]  );
			++n_data;
		}

		std::cout << "data filled correctly" << std::endl;
	}

	void setBackgroundColor( QColor color = QColor(255,255,255)) {
		scene->setBackgroundBrush( QBrush(color) );
	}

	void buildSceneFromData( float w , float h , QColor color = QColor(0,0,255) )
	{
		QPen pen(color);
		QBrush brush(color);
		brush.setStyle(Qt::SolidPattern);
		for( unsigned int r = 0 ; r < repartition.size() ; ++r )
		{
			float x_min = w * (float)(r) / (float)(repartition.size());
			float x_max = w * (float)(r + 1) / (float)(repartition.size());
			float y_max = h * (float)(repartition[r]) / (float)(repartition_max);
			scene->addRect( QRectF(x_min, h-y_max , x_max-x_min , y_max)  ,  pen , brush);
		}
	}

	template< class color_function_t >
	void buildSceneFromDataWithColor( float w , float h )
	{
		for( unsigned int r_it = 0 ; r_it < repartition.size() ; ++r_it )
		{
			float r,g,b;
			float x_center = ((float)(r_it)+0.5f) / (float)(repartition.size());
			color_function_t::getColor(x_center,r,g,b);
			QPen pen(QColor(255*r,255*g,255*b));
			QBrush brush(QColor(255*r,255*g,255*b));
			brush.setStyle(Qt::SolidPattern);
			float x_min = w * (float)(r_it) / (float)(repartition.size());
			float x_max = w * (float)(r_it + 1) / (float)(repartition.size());
			float y_max = h * (float)(repartition[r_it]) / (float)(repartition_max);
			scene->addRect( QRectF(x_min, h-y_max , x_max-x_min , y_max)  ,  pen , brush);
		}
		std::cout << "buildSceneFromDataWithColor done" << std::endl;
	}

	void renderImage( const std::string & filename , int w , int h )
	{
		QImage _image(w,h  ,  QImage::Format_RGB32);

		QPainter painter(&_image);
		scene->render(&painter);

		_image.save(QString::fromStdString(filename));
		std::cout << filename << " saved" << std::endl;
	}
};






template< typename T >
class simpleDataManager{
	std::vector< T > _data;
	T _data_min, _data_max;
	DrawHistogram _histo;
	bool keepMinValue , keepMaxValue;
public:
	simpleDataManager() : _data_min( FLT_MAX ) , _data_max( -FLT_MAX ) , keepMinValue(false) , keepMaxValue(false) { _data.clear(); }
	void add( T * _new_data , int n )
	{
		for( int i = 0 ; i  < n ; ++i )
		{
			_data.push_back(_new_data[i]);
			_data_min = std::min( _data_min , _new_data[i] );
			_data_max = std::max( _data_max , _new_data[i] );
		}
	}
	void add( std::vector< T > const & _new_data )
	{
		for( unsigned int i = 0 ; i  < _new_data.size() ; ++i )
		{
			_data.push_back(_new_data[i]);
			_data_min = std::min( _data_min , _new_data[i] );
			_data_max = std::max( _data_max , _new_data[i] );
		}
	}

	std::vector<T> & getData(){ return _data; }
	T getMinData(){ return _data_min; }
	T getMaxData(){ return _data_max; }

	void forceMinData( T minVal )
	{
		_data_min = minVal;
		keepMinValue = true;
	}

	void renderHistogramImage( const std::string & filename , int w , int h , float relevant , int n_intervalles)
	{
		_histo.clear();
		_histo.fillData( getData() , getMinData() , getMaxData() , relevant , n_intervalles , keepMinValue , keepMaxValue );
		_histo.buildSceneFromData(1,1);
		_histo.renderImage(filename , w , h);
	}
};




template< typename T , class point_t >
class tetraMeshTexture{
	std::vector< T > values;
	std::vector< point_t > vertices;
	std::vector< int > tetras;

	point_t center;
	double scale;
public:
	tetraMeshTexture( const std::vector< point_t > & v , const std::vector< int > & t ) : center(0,0,0) , scale(1.0)
	{
		vertices = v;
		tetras = t;
	}
	void setCenter( const point_t & c ){ center = c; }
	void setScale( double s ){ scale = s; }
	std::vector< T > & getValues(){ return values; }
	void export_vtk( std::string const & filename , std::string const & field_name , std::string const & field_type )
	{
		std::ofstream myfile;
		myfile.open(filename.c_str());
		if (!myfile.is_open())
		{
			std::cout << filename << " cannot be opened" << std::endl;
			return;
		}

		myfile << "# vtk DataFile Version 2.0" << std::endl;
		myfile << filename << std::endl << "ASCII" << std::endl << "DATASET UNSTRUCTURED_GRID" << std::endl << "POINTS " << vertices.size() << " FLOAT" << std::endl;
		for( unsigned int v = 0 ; v < vertices.size() ; ++v )
			myfile << (center + scale*vertices[v]) << std::endl;
		myfile << "CELLS " << tetras.size()/4 << " " << (5 * tetras.size()/4) << std::endl;
		for( unsigned int t = 0 ; t < tetras.size()/4 ; ++t )
			myfile << "4 " << tetras[ 4*t ] << " " << tetras[ 4*t+1 ] << " " << tetras[ 4*t+2 ] << " " << tetras[ 4*t+3 ] << std::endl;
		myfile << "CELL_TYPES " << tetras.size()/4 << std::endl;
		for( unsigned int t = 0 ; t < tetras.size()/4 ; ++t )
			myfile << "10" << std::endl;
	   // myfile << std::endl;
		myfile << "POINT_DATA " << vertices.size() << std::endl;
		myfile << "SCALARS " << field_name << " " << field_type << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
		for( unsigned int i = 0 ; i < vertices.size() ; ++i)
			myfile << values[i] << std::endl;

		myfile.close();
	}

};



template< typename T , class point3d >
class regularGridTexture{
	int nx, ny, nz;
	point3d _origin;
	point3d _spacing;
	std::vector< T > values;
public:
	regularGridTexture( unsigned int x , unsigned int y , unsigned int z ){ nx = x; ny = y; nz = z; }
	void clear(){ values.clear(); }

	void setOrigin( point3d const & o ){ _origin = o; }
	void setSpacing( point3d const & o ){ _spacing = o; }

	void setDim( unsigned int x , unsigned int y , unsigned int z ){ nx = x; ny = y; nz = z; }
	std::vector< T > & getValues( ){ return values; }

	void export_vtk( std::string const & filename , std::string const & field_name , std::string const & field_type , int field_components )
	{
		std::ofstream myfile;
		myfile.open(filename.c_str());
		if (!myfile.is_open())
		{
			std::cout << filename << " cannot be opened" << std::endl;
			return;
		}

		myfile << "# vtk DataFile Version 2.0" << std::endl;
		myfile << filename << std::endl << "ASCII" << std::endl << "DATASET STRUCTURED_POINTS" << std::endl << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
		myfile << "ORIGIN " << _origin << std::endl;
		myfile << "SPACING " << _spacing << std::endl;
		myfile << "POINT_DATA " << nx*ny*nz << std::endl;
		myfile << "SCALARS " << field_name << " " << field_type << " " << field_components << std::endl;
		myfile << "LOOKUP_TABLE default" << std::endl;
		for( unsigned int i = 0 ; i < nx*ny*nz ; ++i)
			myfile << values[i] << std::endl;

		myfile.close();
	}


	void export_vtk( std::string const & filename , std::string const & field_name , std::string const & field_type , int field_components,
					std::vector< std::pair< std::string , double > > const & max_relevant )
	{
		std::vector< T > ordered_values( values.size() );
		for( unsigned int v = 0 ; v < values.size() ; ++v )
			ordered_values[v] = values[v];

		std::sort( ordered_values.begin() , ordered_values.end() );
		for( unsigned int mr = 0 ; mr < max_relevant.size() ; ++mr )
		{
			double max_value = ordered_values[ std::min( (int) ( (double)(values.size()) * max_relevant[mr].second ) , (int)values.size() - 1) ];
			std::string suffix = max_relevant[mr].first;
			std::ofstream myfile;
			myfile.open((filename + suffix + ".vtk").c_str());
			if (!myfile.is_open())
			{
				std::cout << filename << " cannot be opened" << std::endl;
				return;
			}

			myfile << "# vtk DataFile Version 2.0" << std::endl;
			myfile << filename << std::endl << "ASCII" << std::endl << "DATASET STRUCTURED_POINTS" << std::endl << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
			myfile << "ORIGIN " << _origin << std::endl;
			myfile << "SPACING " << _spacing << std::endl;
			myfile << "POINT_DATA " << nx*ny*nz << std::endl;
			myfile << "SCALARS " << field_name << " " << field_type << " " << field_components << std::endl;
			myfile << "LOOKUP_TABLE default" << std::endl;
			for( unsigned int i = 0 ; i < nx*ny*nz ; ++i)
				myfile << std::min(max_value,values[i]) << std::endl;

			myfile.close();
		}
	}
};


#endif // DRAWHISTOGRAM_H
