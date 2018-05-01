
function get_object_top_coord( obj ) {

	var top = obj.offsetTop;
	var height = obj.offsetHeight;
	while( obj=obj.offsetParent ){
		top += obj.offsetTop;
	}

	return top;
}

function get_object_left_coord( obj ) {

	var left=obj.offsetLeft;
	while( obj=obj.offsetParent ){
		left += obj.offsetLeft;
	}

	return left;
}

function mouseCoords(ev){
	if(ev.pageX || ev.pageY){
		return {x:ev.pageX, y:ev.pageY};
	}
	return {
		x:ev.clientX + document.body.scrollLeft - document.body.clientLeft,
		y:ev.clientY + document.body.scrollTop  - document.body.clientTop
	};
}

function M_showTrHint(layer_id, obj, v) {

    var hint_layer = document.getElementById(layer_id);
    var left = obj.offsetLeft;
    var top = get_object_top_coord(obj);
    var obj_height = obj.offsetHeight;

    top += obj_height;
    left -= 10; // because of div.hidable padding-left is 10px

    var style;
    if( style = hint_layer.style ) {

        v=(v=='show')?'visible':(v=='hide')?'hidden':v;

		style.top = top + "px";
		style.left = left + "px";

		style.visibility = v;
    }
}

function M_showTdHint(layer_id, obj, v) {

    var hint_layer = document.getElementById(layer_id);
    var tr_obj = obj.offsetParent;
    var left = tr_obj.offsetLeft;

    var top = get_object_top_coord(obj);
    var obj_height = obj.offsetHeight;

    top += obj_height;
    left -= 10; // because of div.hidable padding-left is 10px

    var style;
    if( style = hint_layer.style ) {

        v=(v=='show')?'visible':(v=='hide')?'hidden':v;

		style.top = top + "px";
		style.left = left + "px";

		style.visibility = v;
    }
}

function showRunningProgramAlert( form_div_id, alert_msg_div_id ) {


    var form_div = document.getElementById(form_div_id);
    var alert_msg_div = document.getElementById(alert_msg_div_id);

    form_div.style.display = 'none';
    alert_msg_div.style.display = 'block';
}

