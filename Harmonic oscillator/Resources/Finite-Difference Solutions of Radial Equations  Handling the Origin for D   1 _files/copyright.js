function addLink() {
    console.log('hey');
    var body_element = document.getElementsByTagName('body')[0];

    var selection;
    selection = window.getSelection();

    // added lines begin
    var htmlDiv = document.createElement("div");
    for (var i = 0; i < selection.rangeCount; ++i) {
        htmlDiv.appendChild(selection.getRangeAt(i).cloneContents()); 
    } 
    var selectionHTML = htmlDiv.innerHTML;    
    // added lines end

    var pagelink = "<br /><br />Read more: http://www.physicsforums.com <br />";
    var copytext = selectionHTML + pagelink; // <------------------- changed line

    var newdiv = document.createElement('div');
    newdiv.style.position='absolute';
    newdiv.style.left='-99999px';
    body_element.appendChild(newdiv);
    newdiv.innerHTML = copytext;
    selection.selectAllChildren(newdiv);
    window.setTimeout(function () {
        body_element.removeChild(newdiv);
    }, 0);
}
document.oncopy = addLink;